
package expressiontable;

import java.io.IOException;
import java.util.HashMap;
import org.apache.commons.math.stat.correlation.SpearmansCorrelation;
import umcg.genetica.io.ExpressionDataset;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author dashazhernakova
 */
public class Coexpression {
    
    public void calculateCoexpression(String fname, String out_fname) throws IOException{
        TextFile out = new TextFile(out_fname, true);
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset(fname);
        double[][] rawData = dataset.getRawData();
        dataset.recalculateHashMaps();
        HashMap<String, Integer> hashProbes = new HashMap(dataset.hashRows);//probes to indices in rawData
        double[] pr1_expr, pr2_expr;
        double cor;
        
        String probe1 = "7_50472431";    
        for (String probe2 : hashProbes.keySet()){
                if (! probe1.equals(probe2)){
                    pr1_expr = rawData[hashProbes.get(probe1)];
                    pr2_expr = rawData[hashProbes.get(probe2)];
                    cor = new SpearmansCorrelation().correlation(pr1_expr, pr2_expr);
                    out.writeln(probe1 + "\t" + probe2 + "\t" + cor);
                }
            }
        //}
	out.close();
    }
    public static void main(String[] args) throws IOException {
        Coexpression c = new Coexpression();
        c.calculateCoexpression("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/tagwise_expression_table_SNP_in_recognition_sequence_tags_excluded.txt.QuantileNormalized.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/coexpression_noPCA_7_50472431");
    }
}
