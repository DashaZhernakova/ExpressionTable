package expressiontable;

import java.io.IOException;
import java.util.HashSet;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author dashazhernakova
 */
public class Joiner {
    DoubleMatrixDataset<String, String> table1;
    DoubleMatrixDataset<String, String> table2;
    public Joiner(String f1, String f2) throws IOException{
        table1 = new DoubleMatrixDataset<String, String>(f1);
        table2 = new DoubleMatrixDataset<String, String>(f2);
    }
    public void addNewProbes(String outFileName) throws IOException{
        TextFile out = new TextFile(outFileName, true);
        HashSet<String> newProbes = new HashSet<String>();
        
        //looking for probes from table2 not present in table1
        for (String probe : table2.rowObjects){
            if (! table1.rowObjects.contains(probe))
                newProbes.add(probe);
        }
        out.close();
    }
    
    public void appendSamples(String outFileName) throws IOException{
        TextFile out = new TextFile(outFileName, true);
        int lineN1 = 0, lineN2 = 0;
        //header
        for (String id : table1.colObjects)
            out.write("\t" + id);
        for (String id : table2.colObjects)
            out.write("\t" + id);
        out.writeln();
        //probes+expression
        for (String probe : table1.rowObjects){
            if (table2.rowObjects.contains(probe)){
                out.write(probe);
                lineN1 = table1.hashRows.get(probe);
                lineN2 = table2.hashRows.get(probe);
                for (int i = 0; i < table1.nrCols;i++){
                    out.write("\t" + table1.rawData[lineN1][i]);
                }
                for (int i = 0; i < table2.nrCols;i++){
                    out.write("\t" + table2.rawData[lineN2][i]);
                }
                out.writeln();
            }
                
        }
        
        out.close();
    }
    
    public void merge(){
        
    }
    public static void main(String[] args) throws IOException {
        /*Joiner j = new Joiner("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Montgomery/expression_table_all.txt.expressedInAllSamples.txt.200genes.sortedByName.txt.QuantileNormalized.Log2Transformed.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Yale+Argonne/expression_table_all_yale+argonne.txt.expressedInAllSamples.txt.genes.txt.QuantileNormalized.Log2Transformed.txt");
        
        j.appendSamples("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Montgomery/Montgomery+Pickrell.genes.txt.QuantileNormalized.Log2Transformed.txt");
         * 
         */
        System.out.println("tfd/sdfgs/sdfgs.txt/fsdf.gz".replaceAll("(\\.gz)?(\\.txt)?$", ""));
    }
}
