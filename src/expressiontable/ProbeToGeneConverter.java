package expressiontable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class ProbeToGeneConverter {
    HashMap<String,String> probe2genes;
    public ProbeToGeneConverter(String annotationFile) throws IOException{
        TextFile annotation = new TextFile(annotationFile, false);
        probe2genes = new HashMap<String, String>();
        String[] els = annotation.readLineElems(TextFile.tab);
        while ((els = annotation.readLineElems(TextFile.tab)) != null)
            probe2genes.put(els[1], els[2]);
        annotation.close();
        
    }
    
    public ProbeToGeneConverter(){}
    
    /*
     * Converts one gene ids X to gene ids Y
     * fname - path to the expression table
     * conversionFname - path to the file of the type X \t Y
     */
    public void convertGeneIdsToGeneNames(String fname, String outFname, String conversionFname) throws IOException{
        HashMap<String, String> conversion = new HashMap<String, String>();
        TextFile conv = new TextFile(conversionFname, false);
        String [] els;
        while ((els = conv.readLineElems(TextFile.tab)) != null){
            conversion.put(els[0], els[1]);
        }
        conv.close();
        
        TextFile table = new TextFile(fname, false);
        TextFile out = new TextFile(outFname, true);
        out.writeln(table.readLine());
        String line;
        int pos = 0, neg = 0;
        while ((line = table.readLine()) != null){
            els = line.split("\t");
            if (conversion.containsKey(els[0])){
                els[0] = conversion.get(els[0]);
                out.writelnTabDelimited(els);
                pos++;
            }
            else
                neg++;
        }
        System.out.println("Successfully converted " + pos + " genes\nNo alternative name found for " + neg + " genes.");
        table.close();
        out.close();
    }
    
    /*
     * writes expression values averaged over all isoforms of a gene
     * fname - expression table
     * unique - write only genes with one isoform
     */
    public void convertProbesToGenesAvg(String fname, String outFname, boolean unique) throws IOException{
        TextFile expr = new TextFile(fname, false);
        
        TextFile out = new TextFile(outFname, true);
        String[] spl;
        String line= expr.readLine(), probe = null, gene = null;
        out.writeln(line);
        TreeMap<String, ArrayList<String>> gene2lines = new TreeMap<String, ArrayList<String>>();
        TreeMap<String, String> gene2avg = new TreeMap<String, String>();
        System.out.println("Converting only single isoform genes? " + unique);
        int numProbes = 0;
        while ((line = expr.readLine()) != null){
            spl = line.split("\t");
            probe = spl[0];
            numProbes ++;
            if (probe2genes.containsKey(probe)){
                gene = probe2genes.get(probe);
             
                ArrayList<String> lines = new ArrayList<String>();
                if (gene2lines.containsKey(gene)) {
                    gene2lines.get(gene).add(line);
                }
                else{
                    lines.add(line);
                    gene2lines.put(gene, lines);
                }
            }
        }
        System.out.println("Overall number of probes processed: " + numProbes);
        System.out.println("Overall number of resulting genes: " + gene2lines.keySet().size());
        
        //Averaging and writing to file
        int size;
        String[] splLine;
        float[] sum;
        String avg; // average gene expression for each sample
        for (Entry <String, ArrayList<String>> e : gene2lines.entrySet()){
            ArrayList<String> lines = e.getValue();
            gene = e.getKey();
            size = lines.size();
            sum = new float[lines.get(0).split("\t").length];
            avg = "";
            if ((size > 1) && (! unique)){ //if more than one isoform for this gene
                for (String s : lines){
                    splLine = s.split("\t");
                    for (int i = 1; i < splLine.length; i++)
                        sum[i]+=Float.parseFloat(splLine[i]); //summing
                }
                out.write(gene);
                for (int i = 1; i < sum.length; i++){ 
                    out.write("\t" + sum[i]/size); //averaging over isoform expression values for current sample
                }
                out.writeln();
                //gene2avg.put(gene, avg);
            }
            else if (size == 1){ //if one isoform
                out.write(gene);//average gene expression = isoform expression
                splLine = lines.get(0).split("\t");
                for (int i = 1; i < splLine.length; i++)
                   out.write("\t" + splLine[i]);
                out.writeln();
            }
        }
        expr.close();
        out.close();
    }
    public static void main(String[] args) throws IOException {
        /*ProbeToGeneConverter c = new ProbeToGeneConverter();
        
        c.convertGeneIdsToGeneNames("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/Pickrell/genes/expression_table.Pickrell.genes.txt.gz.QuantileNormalized.Log2Transformed.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/Pickrell/genes/expression_table.Pickrell.geneNames.txt.gz.QuantileNormalized.Log2Transformed.txt", 
                "/Users/dashazhernakova/Documents/UMCG/hg19/Ids_conversion/Ensembl_v69_geneId2gene.txt");
                */
        ProbeToGeneConverter c = new ProbeToGeneConverter("/Users/dashazhernakova/Documents/UMCG/hg19/annotations/annotation_tag_hg19.txt");
        c.convertProbesToGenesAvg("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/randomSubsets/45samples/expression_table.deepSAGE_tag.45samples.1.txt.gz.QuantileNormalized.Log2Transformed.txt.gz", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/randomSubsets/45samples/expression_table.deepSAGE_tag.45samples.1.txt.gz.QuantileNormalized.Log2Transformed_genes.txt.gz", false);
    }
}
