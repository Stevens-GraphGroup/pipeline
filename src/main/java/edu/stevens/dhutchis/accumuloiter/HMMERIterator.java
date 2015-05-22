package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.data.*;
import org.apache.accumulo.core.iterators.IteratorEnvironment;
import org.apache.accumulo.core.iterators.SortedKeyValueIterator;
import org.apache.commons.lang.SerializationUtils;
import org.apache.hadoop.io.Text;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.RunnableFuture;

/**
 * Calls HMMER native library.
 * Takes in sequences from TseqRaw.
 * Batches into 5000 sequences each.
 * Calls HMMER, which returns an array of bool
 * Packages up results from all 5000 sequences.
 * Emits an entry with all those results.
 *
 * Todo: Timing with Watch
 */
public class HMMERIterator implements SortedKeyValueIterator<Key,Value> {
  private static final Logger log = Logger.getLogger(HMMERIterator.class);

  private SortedKeyValueIterator<Key,Value> source;
  private String hmm_path = "/home/echerin/48.hmm";
  private RangeSet rowRanges = new RangeSet();
  private Iterator<Range> rangeIter;
  private boolean inclusive;
  private Collection<ByteSequence> columnFamilies;
  private Key topKey;
  private Value topValue;
  private int batchSize = 10000;
  private long maxNumBytes = 1000000000;

//  private static final AtomicBoolean loadedNativeLibraries = new AtomicBoolean(false);

  // Load native library
  static {
    // Check standard directories
//    List<File> directories = new ArrayList<File>(Arrays.asList(new File[]{new File("/usr/lib64"), new File("/usr/lib")}));
    // Check in ACCUMULO_HOME location, too
    String envAccumuloHome = System.getenv("ACCUMULO_HOME");
//    if (envAccumuloHome != null) {
//      directories.add(new File(envAccumuloHome + "/lib/native"));
//      directories.add(new File(envAccumuloHome + "/lib/native/map")); // old location, just in case somebody puts it here
//    }
    File libFile = new File(envAccumuloHome+"/lib/ext/Wrap.so");

    // Attempt to load from these directories, using standard names
    if (libFile.exists() && libFile.isFile()) {
      String errMsg = "Tried and failed to load native map library " + libFile;
      try {
        System.load(libFile.getAbsolutePath());
        log.info("Loaded native map shared library " + libFile);
      } catch (Exception | UnsatisfiedLinkError e) {
        log.error(errMsg, e);

        String ldLibraryPath = System.getProperty("java.library.path");
        errMsg = "Tried and failed to load native map library from " + ldLibraryPath;
        try {
          System.loadLibrary("Wrap");
          log.info("Loaded native map shared library from " + ldLibraryPath);
        } catch (Exception | UnsatisfiedLinkError e2) {
          log.error(errMsg, e2);
        }

      }
    } else {
      log.debug("Native map library " + libFile + " not found or is not a file.");
    }
  }


  //  @SuppressWarnings("unchecked")
  private Value hmmerAttachBool(String[] accIDs, String[] rawSeqs) {
    // TODO: at some point later, return the score/probability
    boolean[] booleans = Wrap.seqpass(rawSeqs, hmm_path);

    HashMap<String,String> map = new HashMap<>(rawSeqs.length);
    for (int i = 0; i < rawSeqs.length; i++) {
      map.put(accIDs[i], (booleans[i] ? "1" : "0") /*+ rawSeqs[i]*/);
    }

    byte[] bytes = SerializationUtils.serialize(map);
    return new Value(bytes);
  }

  private static final byte[] SEQ_COL = new Text("seq").copyBytes();

  private void prepareNextEntry() throws IOException {
    topKey = null;
    topValue = null;
    long numBytes=0;
    List<String> accIDs = new ArrayList<>(), rawSeqs = new ArrayList<>();
    Key k = new Key();
    for (int i = 0; numBytes < maxNumBytes && i < batchSize && (source.hasTop() || rangeIter.hasNext()); source.next(), i++) {
      // if finished this range, then goto next range
      while (!source.hasTop() && rangeIter.hasNext())
        source.seek(rangeIter.next(), columnFamilies, inclusive);
      if (!source.hasTop())
        break; // finished ranges; completely done

      // must have "seq" col
      if (!Arrays.equals(source.getTopKey().getColumnQualifierData().getBackingArray(), SEQ_COL))
        continue;
      k.set(source.getTopKey());
      accIDs.add(source.getTopKey().getRow().toString());
      String rawSeq = source.getTopValue().toString();
      numBytes += rawSeq.getBytes().length;
      rawSeqs.add(rawSeq);
    }
    if (rawSeqs.size() > 0) {
//      System.out.printf("hmmerAttachBool: rawSeqs.length= %5d numBytes= %9d Thread= %s\n", rawSeqs.size(), numBytes, Thread.currentThread().getName());
//      long free = Runtime.getRuntime().freeMemory();
//      long max = Runtime.getRuntime().maxMemory();
//      long avail = Runtime.getRuntime().totalMemory();

//      System.out.printf("Runtime: Free %,9d Max %,9d Avail %,d\n",free,(max == Long.MAX_VALUE ? -2 : max),avail);
      topKey = new Key(k);
      topValue = hmmerAttachBool(accIDs.toArray(new String[accIDs.size()]),
               rawSeqs.toArray(new String[rawSeqs.size()]));
    }
  }

  /**
   * Call next() on skvi until getTopKey() advances >= keyToSkipTo (in terms of pk), or until !hasTop().
   * Calls seek() if this takes a while, say greater than 10 next() calls.
   *
   * @return True if advanced to a new key; false if !hasTop().
   */
  static boolean skipUntil(SortedKeyValueIterator<Key, Value> skvi, Key keyToSkipTo,
                           Range seekRange, Collection<ByteSequence> columnFamilies, boolean inclusive) throws IOException {
    /** Call seek() if using this many next() calls does not get us to rowToSkipTo */
    final int MAX_NEXT_ATTEMPT = 10;
    int cnt;
    for (cnt = 0;
         cnt < MAX_NEXT_ATTEMPT && skvi.hasTop() && keyToSkipTo.compareTo(skvi.getTopKey()) > 0;
         cnt++) {
        skvi.next();
    }
    if (skvi.hasTop() && keyToSkipTo.compareTo(skvi.getTopKey()) > 0) {
      // set target range to beginning of pk
      Range skipToRange = new Range(keyToSkipTo, true, null, false)
          .clip(seekRange, true);
      if (skipToRange == null) // row we want to get to does not exist, and it is out of our range
        return false;
      skvi.seek(skipToRange, columnFamilies, inclusive);
    }

    return skvi.hasTop();
  }


  @Override
  public void init(SortedKeyValueIterator<Key, Value> source, Map<String, String> map, IteratorEnvironment iteratorEnvironment) throws IOException {
    this.source = source;
    if (map != null && map.containsKey("hmm_path"))
      hmm_path = map.get("hmm_path");
    if (map != null && map.containsKey("rowRanges"))
      rowRanges.setTargetRanges(GraphuloUtil.d4mRowToRanges(map.get("rowRanges")));
    if (map != null && map.containsKey("batchSize"))
      batchSize = Integer.parseInt(map.get("batchSize"));
    if (map != null && map.containsKey("maxNumBytes"))
      maxNumBytes = Long.parseLong(map.get("maxNumBytes"));
  }

  @Override
  public boolean hasTop() {
    return topKey != null;
  }

  @Override
  public void next() throws IOException {
    prepareNextEntry();
  }

  @Override
  public void seek(Range range, Collection<ByteSequence> columnFamilies, boolean inclusive) throws IOException {
//    this.seekRange = range;
    rangeIter = rowRanges.iteratorWithRangeMask(range);
    this.columnFamilies = columnFamilies;
    this.inclusive = inclusive;
    if (rangeIter.hasNext()) {
    	source.seek(rangeIter.next(), columnFamilies, inclusive);
	    prepareNextEntry();
    }
  }

  @Override
  public Key getTopKey() {
    return topKey;
  }

  @Override
  public Value getTopValue() {
    return topValue;
  }

  @Override
  public SortedKeyValueIterator<Key, Value> deepCopy(IteratorEnvironment iteratorEnvironment) {
    HMMERIterator other = new HMMERIterator();
    other.hmm_path = hmm_path;
    other.source = source.deepCopy(iteratorEnvironment);
    other.rowRanges = rowRanges;
    return other;
  }
}
