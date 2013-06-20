package expressiontable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import umcg.genetica.io.ExpressionDataset;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;


/**
 *
 * @author dashazhernakova
 */
public class Subtable {

    public double calculateAvg(double[] array){
        double avg = 0;
        for (int i = 0 ; i < array.length; i++)
            avg += array[i];
        avg /= array.length;
        return avg;
    }
    public boolean isExpressedInAllSamples(double[] array){
        for (int i = 0 ; i < array.length; i++){
            if (array[i] == 0)
                return false;
        }
        return true;
    }
    public boolean isExpressedInNSamples(double[] array, int N){
        int n = 0;
        for (int i = 0 ; i < array.length; i++){
            if (array[i] > 0)
                n++;
        }
        if (n >= N)
           return true;
        return false;
    }
    public boolean avgExpressionHigherThanThreshold(double[] array, double threshold){
        double avg = 0;
        for (int i = 0 ; i < array.length; i++)
            avg += array[i];
        avg /= array.length;
        if (avg < threshold)
            return false;
        return true;
    }
    
    public void getRandomSubsetOfSamples(String fname, int n, String outFname) throws IOException{
        //reading the expression table
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(fname);
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashSamples = new HashMap(dataset.hashCols);//samples to indices in rawData
        
        List<String> samples = dataset.colObjects;
        HashSet<String> samplesToInclude = new HashSet<String>(n);
        
        Collections.shuffle(samples);
        samplesToInclude.addAll(samples.subList(0, n));
        for (String s : samplesToInclude)
            System.out.println(s);
        dataset = new DoubleMatrixDataset<String, String> (fname, new HashSet(dataset.rowObjects), samplesToInclude);
        
        dataset.recalculateHashMaps();
        dataset.save(outFname);
        
        
    }
    public void getAvgExpression(String fname, String outFname) throws IOException{
        TextFile out = new TextFile(outFname, true);
        
        //reading the expression table
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(fname);
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashProbes = new HashMap(dataset.hashRows);//probes to indices in rawData
        
        //out.writeln("gene\tavg");
        String probe = "";
        double[] line = null;
        int lineNum = 0;
        for ( Entry <String, Integer> e : hashProbes.entrySet()){
            lineNum = e.getValue(); //probe index in rawData
            probe = dataset.rowObjects.get(lineNum);
            line = dataset.getRawData()[lineNum]; //probe expression
            out.writeln(probe + "\t" + calculateAvg(line));
       }
       out.close();
    }
   
    /**
     * Gets top N most expressed probes (N specified by numProbes)
     * @param fname
     * @param outFname
     * @param numProbes
     * @throws IOException 
     */
    public void getMostExpressed(String fname, String outFname, int numProbes) throws IOException{
        TextFile out = new TextFile(outFname, true);
        
        //reading the expression table
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(fname);
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashProbes = new HashMap(dataset.hashRows);//probes to indices in rawData
        
        out.write("\t");
        out.writelnTabDelimited(dataset.colObjects.toArray());
        
        double[] line = null;
        int lineNum = 0;
        HashMap<Integer, Double> probeNumToAvg = new HashMap<Integer, Double>(); //probe indices in rawData to avg expresion
        Sorter s = new Sorter();
        Sorter.ValueComparator bvc =  s.new ValueComparator(probeNumToAvg); //to sort by value (avg expression) rather than by key
        TreeMap<Integer,Double> sorted_probeNumToAvg = new TreeMap<Integer,Double>(bvc); //probeNumToAvg sorted by avg expression

        for ( Entry <String, Integer> e : hashProbes.entrySet()){
            lineNum = e.getValue(); //probe index in rawData
            line = dataset.getRawData()[lineNum]; //probe expression
            probeNumToAvg.put(lineNum, calculateAvg(line));
        }
        sorted_probeNumToAvg.putAll(probeNumToAvg);
        for (Entry<Integer,Double> e : sorted_probeNumToAvg.entrySet()){
            lineNum = e.getKey();
            if (lineNum < numProbes){
                line = dataset.getRawData()[lineNum];
                out.write(dataset.rowObjects.get(lineNum));
                for (int i = 0; i < line.length; i++)
                    out.write("\t" + line[i]);
                out.writeln();
            }
        }
        out.close();
    }
    
    
    public void getExpressedInAllSamples(String fname, String outFname) throws IOException{
        TextFile out = new TextFile(outFname, true);
        
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(fname);
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashProbes = new HashMap(dataset.hashRows);//probes to indices in rawData
        
        out.write("\t");
        out.writelnTabDelimited(dataset.colObjects.toArray());
        
        int counter = 0;
        double[] line = null;
        for ( Entry <String, Integer> e : hashProbes.entrySet()){
            line = dataset.getRawData()[e.getValue()];
            if (isExpressedInAllSamples(line)){
                counter ++;
                out.write(e.getKey());
                for (int i = 0; i < line.length; i++)
                    out.write("\t" + line[i]);
                out.writeln();
            }
        }
        System.out.println("Number of probes expressed in all samples: " + counter);
        out.close();
    }

    /**
     * Writes all probes expressed in at least "percent" % samples
     * @param fname
     * @param outFname
     * @param percent
     * @throws IOException 
     */
    public void getExpressedInNSamples(String fname, String outFname, int percent) throws IOException{
        TextFile out = new TextFile(outFname, true);
        
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(fname);
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashProbes = new HashMap(dataset.hashRows);//probes to indices in rawData
        
        out.write("\t");
        out.writelnTabDelimited(dataset.colObjects.toArray());
        
        int minSamplesExpressed = dataset.nrCols*percent/100;
        int counter = 0;
        double[] line = null;
        for ( Entry <String, Integer> e : hashProbes.entrySet()){
            line = dataset.getRawData()[e.getValue()];
            if (isExpressedInNSamples(line, minSamplesExpressed)){
                counter ++;
                out.write(e.getKey());
                for (int i = 0; i < line.length; i++)
                    out.write("\t" + line[i]);
                out.writeln();
            }
        }
        System.out.println("Number of probes expressed in " + percent + " % of samples: " + counter);
        out.close();
    }
    
     public static void main(String[] args) throws IOException {
         Subtable c = new Subtable();
         /*c.getExpressedInAllSamples("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt", 
                 "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt");
         Sorter s = new Sorter();
         s.sortByAvgExpression("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt", 
                 "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt");
          * 
          */
         //c.getExpressedInAllSamples("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne/yale_argonne_expression_nonnorm_NONZERO.txt", 
        //         "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne/yale_argonne_expression_nonnorm_NONZERO.txt.expressedInAllSamples");
         
         //c.getAvgExpression("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt",
         //        "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/avgExpression.txt");
         c.getRandomSubsetOfSamples("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt", 
                 40,
                 "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt");
                 //"/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/randomSubsets/55samples/expression_table.deepSAGE_tag.55samples.4.txt.gz");
     }
    
  
}
