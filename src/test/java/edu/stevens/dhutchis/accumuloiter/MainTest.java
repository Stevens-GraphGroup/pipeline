/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.client.*;
import org.apache.accumulo.core.client.security.tokens.PasswordToken;
import org.apache.accumulo.core.data.Range;
import org.apache.accumulo.minicluster.MiniAccumuloCluster;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.*;
import java.util.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class MainTest {
    @Rule
    public TemporaryFolder tempFolder = new TemporaryFolder();

    private String username = "root";
    private String password = "secret";
    private static ClientConfiguration myconfig;
    static {
        String instance = "instance";
        String host = "b412srv.ece.stevens-tech.edu:2181";
        int timeout = 10000; // 10 seconds
        myconfig = ClientConfiguration.loadDefault().withInstance(instance).withZkHosts(host).withZkTimeout(timeout);
    }

    public Connector connectToAccumulo() throws AccumuloSecurityException, AccumuloException {
        Instance instance = new ZooKeeperInstance(myconfig);
        return instance.getConnector(username, new PasswordToken(password));
    }


    private void printList(Collection<?> list, String prefix) {
        System.out.print(prefix+": ");
        for (Object o : list) {
            System.out.print(o+", ");
        }
        System.out.println();
    }

    public MainTest() {

    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    //@Test
    public void testMini() throws Exception {

        File tempDir = tempFolder.newFolder();
        MiniAccumuloCluster accumulo = new MiniAccumuloCluster(tempDir, "password");
        accumulo.start(); // doesn't work on Dylan's computer for some reason.  The OS closes the Zookeeper connection.

        Instance instance = new ZooKeeperInstance(accumulo.getInstanceName(), accumulo.getZooKeepers());
        Connector conn = instance.getConnector("root", new PasswordToken("password"));

        //  innerTest(instance, conn);

        accumulo.stop();
        tempDir.delete();
    }

    @Test
    public void testNormal() throws Exception {
        Connector conn = connectToAccumulo();
        innerTest(conn);
    }


    @Test
    public void lightTest() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();
        PrintWriter writer = new PrintWriter("lightTest" + new Date( ).toString() + ".txt", "UTF-8");
        String tInput = "taxonomy|Bacteria"; //; Firmicutes; Clostridia"; //contains about 500,000 seqs
        System.out.println(wrap.taxToRaw(conn, tInput, 40000, 40000, 20000000, 3));
    }


    @Test
    public void sameBatchOneThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs

        System.out.println("Same batch size, One Thread");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");


        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, iterBatchSize, iterBatchSize, 1000000000,1));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    @Test
    public void sameBatchTwoThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs

        System.out.println("Same batch size, Two Threads");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");


        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, iterBatchSize, iterBatchSize, 1000000000,2));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    @Test
    public void sameBatchThreeThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs

        System.out.println("Same batch size, Three Threads");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");

        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, iterBatchSize, iterBatchSize, 1000000000,3));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    @Test
    public void twiceBatchOneThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs
        System.out.println("Accession batch size twice iterator batch size, One Threads");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");

        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(2*iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, 2*iterBatchSize, iterBatchSize, 1000000000,1));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    @Test
    public void twiceBatchTwoThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs
        System.out.println("Accession batch size twice iterator batch size, Two Threads");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");
        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(2*iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, 2*iterBatchSize, iterBatchSize, 1000000000,2));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    @Test
    public void twiceBatchThreeThread() throws Exception
    {
        Connector conn = connectToAccumulo();

        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";
        Wrap wrap = new Wrap();

        //String tInput = "taxonomy|Bacteria; Firmicutes; Clostridia"; //contains about 500,000 seqs
        String tInput = "taxonomy|Bacteria; Proteobacteria; Epsilonproteobacteria"; //contains about 250,000 seqs
        System.out.println("Accession batch size twice iterator batch size, Three Threads");
        System.out.println("acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time");

        for(int iterBatchSize = 10000; iterBatchSize <= 100000; iterBatchSize+=10000)
        {
            System.out.println(Integer.toString(2*iterBatchSize) + " " + Integer.toString(iterBatchSize) + " " + wrap.taxToRaw(conn, tInput, 2*iterBatchSize, iterBatchSize, 1000000000,3));
            //acc ids batch size; iterator batch size; how many of taxon was in database; total scan time; total compute time
        }
    }

    private void innerTest( Connector conn) throws Exception {
        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";

        Wrap wrap = new Wrap();
        File f = new File("/home/echerin/ppp/pipeline/taxons.txt");
        java.util.Scanner s = new Scanner(f);
        List<String> used = new ArrayList<>();

        //PrintWriter writer = new PrintWriter("graph" + new Date( ).toString() + ".txt", "UTF-8");

        while(s.hasNextLine()) //assume taxonomy| ... has no spaces in family/genus names
        {

            taxon = s.nextLine();
            //System.out.println("We are in outer while loop and using: " + taxon);

            java.util.Scanner wScan = new Scanner(taxon);
            String tInput = "taxonomy|";

            boolean first = true;
            while(wScan.hasNext())
            {
                if(first)
                {
                    tInput += wScan.next();
                    first = false;
                }
                else
                {
                    tInput += "; " + wScan.next();
                }

                //System.out.println("scanning with : " + tInput);
                if(!used.contains(tInput))
                {
                    System.out.println(wrap.taxToRawCompute(conn, tInput,50000,50000, 200000000,3));
                    used.add(tInput);
                }
            }

        }

        //writer.close();

    }


}
