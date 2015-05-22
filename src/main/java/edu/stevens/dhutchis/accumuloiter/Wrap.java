
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.client.*;
import org.apache.accumulo.core.client.Scanner;
import org.apache.accumulo.core.data.Key;
import org.apache.accumulo.core.data.Mutation;
import org.apache.accumulo.core.data.Value;
import org.apache.accumulo.core.iterators.Combiner;
import org.apache.accumulo.core.iterators.IteratorUtil;
import org.apache.accumulo.core.iterators.LongCombiner;
import org.apache.accumulo.core.iterators.user.SummingCombiner;
import org.apache.accumulo.core.security.Authorizations;
import org.apache.commons.lang.SerializationUtils;
import org.apache.hadoop.io.Text;
import org.junit.Assert;
import org.apache.accumulo.core.data.Range;

import java.io.PrintWriter;
import java.util.*;
import java.util.Map.Entry;

/**
 * @author dhutchis
 */
public class Wrap {

  public synchronized static native boolean[] seqpass(String[] n, String hmm_path);

  private final String tableName = "TestTableIterator";
  private final String columnFamily = "";
  private final String columnVisibility = "";

  private void printList(Collection<?> list, String prefix) {
    System.out.print(prefix + ": ");
    for (Object o : list) {
      System.out.print(o + ", ");
    }
    System.out.println();
  }

  public static void main(String[] args) {
    System.out.println("HI im in main i updated again again");
  }

  void xin(Connector conn) throws AccumuloSecurityException, AccumuloException, TableNotFoundException {

    System.out.println("entasdfax to acc");
    // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
    String TseqRaw = "TseqT";
    int numThreads = 1;
    Scanner scan = conn.createScanner(TseqRaw, Authorizations.EMPTY);
    scan.setRange(new Range("taxonomy", "taxonomy~"));
    Set<String> set = new HashSet<String>();

    for (Map.Entry<Key, Value> entry : scan) {

      String en = entry.getKey().toString();

      if (en.contains(";")) {
        set.add(en.substring(9, en.indexOf(';')));
      }

    }
    scan.close();

    for (String s : set) {
      System.out.println(s);
    }

  }

  public List<Range> taxToAcc(Connector conn, List<String> taxaList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException {
    System.out.println("entered tax to acc");
    // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
    String TseqRaw = "TseqT";
    Scanner scan = conn.createScanner(TseqRaw, Authorizations.EMPTY);
    scan.setRange(new Range("taxonomy|Bacteria; Cyanobacteria", "taxonomy|Bacteria; Cyanobacteria~"));

    List<Range> accList = new ArrayList<>();

    // Do the scan

    for (Map.Entry<Key, Value> entry : scan) {
      String acc = entry.getKey().getColumnQualifier().toString();
      // System.out.println("row is: " +entry.getKey().getRow().toString()+" acc is : " + acc );
      accList.add(new Range(acc));
    }
    scan.close();

    return accList;
  }

  public String taxToRaw(Connector conn, String taxon, int accIdSize, int iterBatchSize, long iterGPUMaxNumBytes, int threadNumber) throws AccumuloSecurityException, AccumuloException, TableNotFoundException {
    long startTime = System.currentTimeMillis();

    final int batchSize = accIdSize;// increase size
    //System.out.println("entered tax to raw");
    // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
    final String TseqT = "TseqT";
    final String TseqRaw = "TseqRaw";
    int numThreads = threadNumber;

    Scanner scan = conn.createScanner(TseqT, Authorizations.EMPTY);
    scan.setRange(new Range(taxon, taxon + "~"));

    BatchScanner batScan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);
    batScan.setRanges(Collections.singleton(new Range()));
    final String hmm_path = "/home/echerin/2405.hmm";

    List<Range> accList = new ArrayList<>();
    Map<String, String> rawSeqMap = new HashMap<String, String>();

    long computeTime = 0; //holds total time it takes to filter the raw sequences
    long startComputeTime;

    // Scan TseqT
    int counter = 0;
    for (Map.Entry<Key, Value> entry : scan) {
      String acc = entry.getKey().getColumnQualifier().toString();
      accList.add(new Range(acc));
      counter++;
      if (counter % batchSize == 0)// when we have enough ranges
      {
        // BatchScan TseqRaw with accession IDs from TseqT.
        startComputeTime = System.currentTimeMillis();
        scanTseqRaw(batScan, hmm_path, accList,iterBatchSize, iterGPUMaxNumBytes);
        computeTime += System.currentTimeMillis() - startComputeTime;
        accList.clear();
      }
    }
    scan.close();
    //System.out.println("entering last batch");
    if (!accList.isEmpty()) {
      startComputeTime = System.currentTimeMillis();
      scanTseqRaw(batScan, hmm_path, accList,iterBatchSize, iterGPUMaxNumBytes);
      computeTime += System.currentTimeMillis() - startComputeTime;
      accList.clear();
    }
    batScan.close();

    long totalTime = System.currentTimeMillis() - startTime;
    return Integer.toString(counter) + " " + Long.toString(totalTime - computeTime) + " " + Long.toString(computeTime)+ " " + taxon;
    //taxon string;  how many of taxon was in database; total scan time; total compute time
  }

  //taken from last commit on April 15th
  public String taxToRawCompute(Connector conn, String taxon, int accIdSize, int iterBatchSize, long iterGPUMaxNumBytes, int threadNumber) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
  {
    long startTime = System.currentTimeMillis();

    System.load("/home/echerin/ppp/pipeline/src/main/java/edu/stevens/dhutchis/accumuloiter/Wrap.so");
    String hmm_path = "/home/echerin/48.hmm";


    int batchSize = accIdSize;
    //System.out.println("entered tax to raw");
    // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
    String TseqT = "TseqT";
    String TseqRaw = "TseqRaw";
    int numThreads = threadNumber;

    Scanner scan = conn.createScanner(TseqT, Authorizations.EMPTY);
    scan.setRange(new Range(taxon, taxon + "~"));

    BatchScanner batScan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);

    List<Range> accList = new ArrayList<>();
    Map<String,String> rawSeq = new HashMap<String,String>();

    long computeTime = 0; //holds total time it takes to filter the raw sequences
    long startComputeTime;

    // Do the scan
    int counter = 0;


    for(Map.Entry<Key,Value> entry : scan)
    {
      if(counter%batchSize == 0 && counter !=0)// when we have enough ranges
      {
        long numBytes=0;
        batScan.setRanges(accList);

        for(Map.Entry<Key,Value> batEntry : batScan)
        {
          String seq = batEntry.getValue().toString();
          String mykey = batEntry.getKey().toString();
          rawSeq.put(mykey, seq);

          numBytes += seq.getBytes().length;
          if(numBytes > iterGPUMaxNumBytes)// break up the batch if the
          {
            startComputeTime = System.currentTimeMillis();
            Wrap.seqpass( rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
            computeTime += System.currentTimeMillis() - startComputeTime;
            rawSeq.clear();
            numBytes = 0;
          }
        }
        if(!rawSeq.isEmpty())
        {
          startComputeTime = System.currentTimeMillis();
          Wrap.seqpass(rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
          computeTime += System.currentTimeMillis() - startComputeTime;
          rawSeq.clear();
        }
        accList.clear();
      }

      String acc = entry.getKey().getColumnQualifier().toString();

      accList.add(new Range(acc));

      counter++;
    }


    //scan the last batch
    if(!accList.isEmpty())
    {

      long numBytes=0;
      batScan.setRanges(accList);

      for(Map.Entry<Key,Value> batEntry : batScan)
      {
        String seq = batEntry.getValue().toString();
        String mykey = batEntry.getKey().toString();
        rawSeq.put(mykey, seq);

        numBytes += seq.getBytes().length;
        if(numBytes > iterGPUMaxNumBytes)// break up the batch if the
        {
          startComputeTime = System.currentTimeMillis();
          Wrap.seqpass( rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
          computeTime += System.currentTimeMillis() - startComputeTime;
          rawSeq.clear();
          numBytes = 0;
        }
      }
      if(!rawSeq.isEmpty())
      {
        startComputeTime = System.currentTimeMillis();
        Wrap.seqpass(rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
        computeTime += System.currentTimeMillis() - startComputeTime;
        rawSeq.clear();
      }
      accList.clear();

    }

    batScan.close();
    scan.close();

    long totalTime = System.currentTimeMillis() - startTime;

    return Integer.toString(counter) + " " + Long.toString(totalTime - computeTime) + " " + Long.toString(computeTime)+ " " + taxon;
  }

  @SuppressWarnings("unchecked")
  private static void scanTseqRaw(BatchScanner batScan, String hmm_path, Collection<Range> accList, int iterBatchSize, long iterGPUMaxNumBytes) {
    //                batScan.setRanges(accList);
    batScan.clearScanIterators();
    Map<String, String> options = new HashMap<>();
    options.put("hmm_path", hmm_path);
    options.put("rowRanges", GraphuloUtil.rangesToD4MString(accList));
    options.put("batchSize",Integer.toString(9999999));
    options.put("maxNumBytes", Long.toString(iterGPUMaxNumBytes));
    IteratorSetting itset = new IteratorSetting(18, HMMERIterator.class, options);
    batScan.addScanIterator(itset);

    for (Map.Entry<Key, Value> batEntry : batScan) {
//      System.out.println("A Entry: " + batEntry.getKey() + " -> " + batEntry.getValue());
//      HashMap<String, String> map1 = (HashMap<String, String>) SerializationUtils.deserialize(batEntry.getValue().get());

//      for (Map.Entry<String, String> accToEncodedRawSeq : map1.entrySet()) {
////        String accID = accToEncodedRawSeq.getKey();
////        String tmp = accToEncodedRawSeq.getValue();
////        boolean b = tmp.charAt(0) != '0';
////        String rawSeq = tmp.substring(1);
////        writer.append("accID=" + accID + "  b=" + b + "  rawSeq=" + rawSeq.length() + " chars");// do something with accID, b, rawSeq
////                        System.out.println("accID="+accID+"  b="+b+"  rawSeq="+rawSeq.length()+" chars");
//      }
    }
  }


  /**
   * input is a list of Accession numbers ///used to be accession now its just con
   * /   output is a map from accession numbers to sequences
   */
  public Map<String, String> accToRaw(Connector conn, List<Range> accessionList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException {

    System.out.println("entered acc to raw");
    // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
    String TseqRaw = "TseqRaw";
    int numThreads = 1;
    BatchScanner batScan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);

        /*
        List<Range> accessionRanges = new ArrayList<>(accessionList.size());
        for (String accession : accessionList)
        {
            accessionRanges.add(new Range(accession));
        }
        */


    Map<String, String> rawSeq = new HashMap<String, String>(accessionList.size());
    // Do the scan
    System.out.println("I am about to scan");

    batScan.setRanges(accessionList);
    for (Map.Entry<Key, Value> entry : batScan) {
      String seq = entry.getValue().toString();
      String mykey = entry.getKey().toString();
      //System.out.println("seq is: "  + seq + " key is : " + mykey);
      rawSeq.put(mykey, seq);
    }
    batScan.close();
    //rawSeq.values().toArray();
    return rawSeq;
  }


  public void testIter(Connector conn) throws AccumuloException, AccumuloSecurityException, TableNotFoundException, TableExistsException {

    System.out.println("I am in testIter");

    // check results
    Scanner scan = conn.createScanner("TseqDegT", Authorizations.EMPTY);
    scan.setRange(new Range("|Bacteria; Cyanobacteria"));
    // System.out.println("Scanner range: "+scan.getRange());
    int i = 0;
    for (Entry<Key, Value> entry : scan) {
      System.out.println(entry);
      // Assert.assertEquals("12", entry.getValue().toString());
      i++;
    }
    System.out.println("Number of results is: " + i);
    // try a custom iterator...
  }

}
