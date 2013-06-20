package expressiontable;

import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author dashazhernakova
 */
public class Sorter {
    public double calculateAvg(double[] array){
        double avg = 0;
        for (int i = 0 ; i < array.length; i++)
            avg += array[i];
        avg /= array.length;
        return avg;
    }
    
    public void sortByAvgExpression(String fname, String outFname) throws IOException{
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
        ValueComparator bvc =  new ValueComparator(probeNumToAvg); //to sort by value (avg expression) rather than by key
        TreeMap<Integer,Double> sorted_probeNumToAvg = new TreeMap<Integer,Double>(bvc); //probeNumToAvg sorted by avg expression

        for ( Entry <String, Integer> e : hashProbes.entrySet()){
            lineNum = e.getValue(); //probe index in rawData
            line = dataset.getRawData()[lineNum]; //probe expression
            probeNumToAvg.put(lineNum, calculateAvg(line));
        }
        sorted_probeNumToAvg.putAll(probeNumToAvg);
        for (Entry<Integer,Double> e : sorted_probeNumToAvg.entrySet()){
            lineNum = e.getKey();
            line = dataset.getRawData()[lineNum];
            out.write(dataset.rowObjects.get(lineNum));
            for (int i = 0; i < line.length; i++)
                out.write("\t" + line[i]);
            out.writeln();
        }
        out.close();
    }
    
    public void sortByProbeName(String fname, String outFname) throws IOException{
        TextFile in = new TextFile(fname, false);
        TextFile out = new TextFile(outFname, true);
        
        String line = in.readLine(), probe;
        out.write("\t");
        out.writeln(line);
        
        TreeMap <String, String> probe2expr = new TreeMap<String, String>();
        
        while ( (line = in.readLine()) != null){
            probe = line.split("\t")[0];
            probe2expr.put(probe, line);
        }
        
        for (String pr : probe2expr.keySet()){
            out.writeln(pr + "\t" + probe2expr.get(pr));
        }
        in.close();
        out.close();
    }
    
    public  class ValueComparator implements Comparator<Integer> {

        Map<Integer, Double> base;
        public ValueComparator(Map<Integer, Double> base) {
            this.base = base;
        }


        @Override
        public int compare(Integer a, Integer b) {
            return base.get(b).compareTo(base.get(a));
        }
    }
    public static void main(String[] args) throws IOException {
        Sorter s = new Sorter();
         s.sortByAvgExpression("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Sebo/expression_table_normByGeneLength.txt", 
                 "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Sebo/expression_table_normByGeneLength_sorted.txt");
    }
}

