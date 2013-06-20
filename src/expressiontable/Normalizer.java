package expressiontable;

import java.io.IOException;
import umcg.genetica.io.ExpressionDataset;
import umcg.genetica.math.matrix.DoubleMatrixDataset;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Log2Transform;
import umcg.genetica.math.stats.QuantileNormalization;

/**
 *
 * @author dashazhernakova
 */
public class Normalizer {
    public void normalize(String expressionFile) throws IOException{
        DoubleMatrixDataset<String, String> dataset = new DoubleMatrixDataset<String, String>(expressionFile);
	double[][] rawData = dataset.getRawData();
	String fileNamePrefix = expressionFile;
	
        
        QuantileNormalization.quantilenormalize(rawData);
//        
	DoubleMatrixDataset<String, String> datasetNormalized = new DoubleMatrixDataset<String, String> (dataset.nrRows, dataset.nrCols);

	datasetNormalized.rowObjects = dataset.rowObjects;
	datasetNormalized.colObjects = dataset.colObjects;
	datasetNormalized.setRawData(rawData);
	fileNamePrefix += ".QuantileNormalized";
	datasetNormalized.save(fileNamePrefix + ".txt.gz");
	datasetNormalized = null;


	Log2Transform.log2transform(rawData);

	datasetNormalized = new DoubleMatrixDataset<String, String>(dataset.nrRows, dataset.nrCols);
	datasetNormalized.rowObjects = dataset.rowObjects;
	datasetNormalized.colObjects = dataset.colObjects;
	datasetNormalized.setRawData(rawData);
	fileNamePrefix += ".Log2Transformed";
	datasetNormalized.save(fileNamePrefix + ".txt.gz");
	datasetNormalized = null;

	System.out.println("Standardizing probe mean and standard deviation");
	for (int p = 0; p < dataset.nrRows; p++) {
	    double mean = Descriptives.mean(rawData[p]);
	    double stdev = Math.sqrt(Descriptives.variance(rawData[p], mean));
	    for (int s = 0; s < dataset.nrCols; s++) {
		rawData[p][s] -= mean;
	    }
	}

        dataset.setRawData(rawData);
	fileNamePrefix += ".ProbesCentered";
	dataset.save(fileNamePrefix + ".txt.gz");

	System.out.println("- Standardizing sample mean and standard deviation");
	for (int s = 0; s < dataset.nrCols; s++) {
	    double[] vals = new double[dataset.nrRows];
	    for (int p = 0; p < dataset.nrRows; p++) {
		vals[p] = dataset.getRawData()[p][s];
	    }
	    double mean = Descriptives.mean(vals);
	    for (int p = 0; p < dataset.nrRows; p++) {
		vals[p] -= mean;
	    }
	    double var = Descriptives.variance(vals, mean);
	    double stdev = Math.sqrt(var);
	    for (int p = 0; p < dataset.nrRows; p++) {
		dataset.getRawData()[p][s] = (vals[p] / stdev);
	    }
	}

	datasetNormalized = new DoubleMatrixDataset<String, String>(dataset.nrRows, dataset.nrCols);
	datasetNormalized.rowObjects = dataset.rowObjects;
	datasetNormalized.colObjects = dataset.colObjects;
	datasetNormalized.setRawData(rawData);
	fileNamePrefix += ".SamplesZTransformed";
	datasetNormalized.save(fileNamePrefix + ".txt.gz");
	datasetNormalized = null;
    }
    public static void main(String[] args) throws IOException {
        Normalizer n = new Normalizer();
        n.normalize("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Yale+Argonne/expression_table_all_yale+argonne.txt.expressedInAllSamples.txt.genes.txt");
    }
}
