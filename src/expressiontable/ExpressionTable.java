/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package expressiontable;

import java.io.IOException;
//import eqtlmappingpipeline.normalization.Normalizer;
/**
 *
 * @author dashazhernakova
 */
public class ExpressionTable {

    /**
     * @param args the command line arguments
     */
    public static void usage(){
        System.out.println("--mode\n\t"
                + "ProbeToGeneConverter\n\t"
                + "getExpressedInAllSamples\n\t"
                + "getTopExpressed\n\t"
                + "sort\n\t"
                + "normalize");
    }
    public static void main(String[] args) throws IOException {
        String lincRNA = "/Users/dashazhernakova/Documents/UMCG/lincRNA/annotation_lincRNA_hg19_toGenes.txt",
                transcr = "/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt";
        Subtable sub = new Subtable();
        Sorter sorter = new Sorter();
        Normalizer norm = new Normalizer();
        //Normalizer norm = new Normalizer();
        
        String arg, val, in = null, mode = null, out = null;
        
        int i = 0;
        for (i = 0; i < args.length; i++) {
	    arg = args[i];
	    val = null;

	    if (i + 1 < args.length) {
		val = args[i + 1];
	    }

	    if (arg.equals("--mode")) {
		mode = val;
                //System.out.println("mode");
                break;
	    }
            
	}
         if (mode == null) {
	    System.out.println("ERROR: Please supply --mode");
            usage();
         }
        else if (mode.equals("ProbeToGeneConverter")){
            String annot = null;
            boolean unique = false;
            for (int j = i; j < args.length; j++){ 
                arg = args[j];
                val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) 
                    in = val;
                if (arg.equals("--annot")) {
                    annot = val;
                }
                if (arg.equals("--out")) 
                    out = val;
                if (arg.equals("--unique"))
                    unique = Boolean.valueOf(val);
            }
            if (out == null)
                out = in.replaceAll("(\\.gz)?(\\.txt)?$", "") + ".genes.txt";
            if ( (in == null ) || (annot == null ))
                 System.out.println("Not enough arguments!!!");
            System.out.println("Converting to gene ids: \n\texpression table " + in + "\n\tunique " + unique + "\n\tannotation " + annot);
            ProbeToGeneConverter converter = new ProbeToGeneConverter(annot);
            converter.convertProbesToGenesAvg(in, out, unique);
        }
         else if (mode.equals("getExpressedInAllSamples")){
             
             for (int j = i; j < args.length; j++){ 
                arg = args[j];
                val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) 
                    in = val;
                if (arg.equals("--out")) 
                    out = val;
                
            }
            if (out == null)
                out = in.replaceAll("(\\.gz)?(\\.txt)?$", "") + ".expressedInAllSamples.txt";
            System.out.println("\nGetting probes expressed in all samples from " + in);
            sub.getExpressedInAllSamples(in, out);
        }
        else if (mode.equals("getTopExpressed")){
            int n = 0;
            for (int j = i; j < args.length; j++){ 
                arg = args[j];
                val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) 
                    in = val;
                if (arg.equals("--out")) 
                    out = val;
                if (arg.equals("--n")) 
                    n = Integer.parseInt(val);
            }
            if (out == null)
                out = in.replaceAll("(\\.gz)?(\\.txt)?$", "") + ".top" + n;
            System.out.println("Getting top " + n + " expressed genes/transcripts... from " + in);
            sub.getMostExpressed(in, out, n);
        }
        
        else if (mode.equals("sort")){
            String by = null;
            for (int j = i; j < args.length; j++){ 
                arg = args[j];
                val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) 
                    in = val;
                if (arg.equals("--out")) 
                    out = val;
                if (arg.equals("--by")) 
                    by = val;
            }
            if (out == null)
                out = in + ".sorted";
            if (by.equals("name")){
                System.out.println("Sorting " + in + " by probe name...");
                if (out == null)
                    out = in.replaceAll("(\\.gz)?(\\.txt)?$", "") + ".sortedByName";
                sorter.sortByProbeName(in, out);
        
            }
            else if (by.equals("expression")){
                System.out.println("Sorting " + in + " by average expression...");
                if (out == null)
                    out = in.replaceAll("(\\.gz)?(\\.txt)?$", "") + ".sortedByExpr";
                sorter.sortByAvgExpression(in, out);
            }
            else
                 System.out.println("Wrong \"by\" parameter");
        }
        else if (mode.equals("normalize")){
            
            for (int j = i; j < args.length; j++){ 
                arg = args[j];
                val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) 
                    in = val;
             }
            System.out.println("Normalizing " + in);
            norm.normalize(in);
        }
        else{
             System.out.println("Wrong mode!");
             usage();
        }
        
    }
}
