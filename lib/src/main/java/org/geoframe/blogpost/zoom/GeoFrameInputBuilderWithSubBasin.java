package org.geoframe.blogpost.zoom;

import static org.hortonmachine.hmachine.modules.network.netnumbering.OmsNetNumbering.OMSNETNUMBERING_desiredAreaDelta_DESCRIPTION;
import static org.hortonmachine.hmachine.modules.network.netnumbering.OmsNetNumbering.OMSNETNUMBERING_desiredArea_DESCRIPTION;
import static org.hortonmachine.hmachine.modules.network.netnumbering.OmsNetNumbering.OMSNETNUMBERING_inPoints_DESCRIPTION;

import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import javax.media.jai.iterator.RandomIterFactory;
import javax.media.jai.iterator.WritableRandomIter;

import org.apache.commons.io.FileUtils;
import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.data.DataUtilities;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.geometry.jts.ReferencedEnvelope;
import org.geotools.process.vector.UnionFeatureCollection;
import org.hortonmachine.gears.libs.exceptions.ModelsRuntimeException;
import org.hortonmachine.gears.libs.modules.HMConstants;
import org.hortonmachine.gears.libs.modules.HMModel;
import org.hortonmachine.gears.utils.RegionMap;
import org.hortonmachine.gears.utils.coverage.CoverageUtilities;
import org.hortonmachine.gears.utils.features.FeatureUtilities;
import org.hortonmachine.gears.utils.files.FileUtilities;
import org.hortonmachine.hmachine.modules.basin.rescaleddistance.OmsRescaledDistance;
import org.hortonmachine.hmachine.modules.demmanipulation.markoutlets.OmsMarkoutlets;
import org.hortonmachine.hmachine.modules.network.netnumbering.OmsNetNumbering;
import org.hortonmachine.modules.GeoframeInputsBuilder;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import oms3.annotations.Author;
import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.Label;
import oms3.annotations.License;
import oms3.annotations.Name;
import oms3.annotations.Status;
import oms3.annotations.UI;
import oms3.annotations.Unit;

/**
 * TODO out into @Description annotation.
 * 
 * 
 * 
 * 
 * 
 * @author Daniele Andreis
 *
 */

@Description("Module to prepare input data for the Geoframe modelling environment.")
@Author(name = "Daniele Andreis", contact = "daniele.andreis@unitn.it")
@Keywords("geoframe")
@Label(HMConstants.HYDROGEOMORPHOLOGY)
@Name("_GeoframeInputsBuilderWithSubBasin")
@Status(5)
@License("General Public License Version 3 (GPLv3)")
public class GeoFrameInputBuilderWithSubBasin extends HMModel {
	private static final String BASIN_CSV_SEPARATOR = ";";
	private static final String TOPOLOGY_OUTLET = " 0";
	private static final String WHITE_SPACE = " ";
	private static final String BASINID_KEY = "basinid";
	private static final String TOPOLOGY_CSV_PATH = "topology.csv";
	private static final String NETWORK_SHP_FILE_PATH = "network_complete.shp";
	private static final String SUBBASINS_COMPLETE_SHP_FILE_PATH = "subbasins_complete.shp";
	private static final String SUBBASINS_CSV_PATH = "subbasins.csv";
	private static final String NEW_LINE = "\n";
	@Description("Id basin to refine")
	@In
	public boolean pGeoFrameAlreadyDone = false;
	@Description("Id basin to refine")
	@In
	public int[] pBasinToRefineIds = null;
	@Description(OMSNETNUMBERING_desiredArea_DESCRIPTION)
	@Unit("m2")
	@In
	public Double pDesiredArea = null;

	@Description(OMSNETNUMBERING_desiredAreaDelta_DESCRIPTION)
	@Unit("%")
	@In
	public Double pDesiredAreaDelta = null;

	@Description("Input pitfiller raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inPitfiller = null;

	@Description("Input flowdirections raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inDrain = null;

	@Description("Input tca raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inTca = null;

	@Description("Input network raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inNet = null;

	@Description("Input skyview factor raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inSkyview = null;

	@Description("Input numbered basins raster map.")
	@UI(HMConstants.FILEIN_UI_HINT_RASTER)
	@In
	public String inBasins = null;

	@Description("Optional input lakes vector map.")
	@UI(HMConstants.FILEIN_UI_HINT_VECTOR)
	@In
	public String inLakes = null;

	@Description("The geoframe topology file, mandatory in case of lakes.")
	@UI(HMConstants.FILEIN_UI_HINT_GENERIC)
	@In
	public String inGeoframeTopology = null;

	@Description(OmsRescaledDistance.OMSRESCALEDDISTANCE_pRatio_DESCRIPTION)
	@In
	public double pRatio = 50;

	@Description("Output folder for the geoframe data preparation")
	@UI(HMConstants.FOLDERIN_UI_HINT)
	@In
	public String outFolder = null;
	@Description("Output folder for the geoframe data preparation")
	@UI(HMConstants.FOLDERIN_UI_HINT)
	@In
	public String outFolderWithSmallBasin = null;
	@Description("Output folder for the geoframe data preparation")
	@UI(HMConstants.FOLDERIN_UI_HINT)
	@In
	public String newFolder = null;

	@Description(OMSNETNUMBERING_inPoints_DESCRIPTION)
	@UI(HMConstants.FILEIN_UI_HINT_VECTOR)
	@In
	public String inPoints = null;
	@Description("Run on each subbasin")
	@UI(HMConstants.FILEIN_UI_HINT_GENERIC)
	@In
	public Boolean inRunAll = false;

	/**
	 * overwrite files
	 */
	private boolean doOverWrite = true;

	private String outMarkOutlet = null;

	/**
	 * String used to build the path for zoomed basin.
	 */
	private static String basinPathSeparator = "000";
	/**
	 * the working directory (basins to zoom)
	 */
	private String subBasinPath;
	/**
	 * Collection used to merge the network in normal basins and network in zoomed
	 * basins,
	 */
	private SimpleFeatureCollection newNetworkFC;
	/**
	 * Collection used to merge the basins boundary(normal plus zoomed basins)
	 */
	private SimpleFeatureCollection newBasinFC;
	/**
	 * Basins id list,
	 */
	private List<String> bigBasinList;

	@Execute
	public void process() throws Exception {
		// if not already done run on whole region.
		if (!this.pGeoFrameAlreadyDone) {
			this.prepareBigBasin();
		}
		// read the shape file of the whole region
		String netBasin = Paths.get(outFolder, NETWORK_SHP_FILE_PATH).toString();
		newNetworkFC = getVector(netBasin);
		String subBasin = Paths.get(outFolder, SUBBASINS_COMPLETE_SHP_FILE_PATH).toString();
		newBasinFC = getVector(subBasin);
		String bigBasinCsv = Paths.get(outFolder, SUBBASINS_CSV_PATH).toString();
		bigBasinList = FileUtilities.readFileToLinesList(bigBasinCsv);
		// iterate over the selected basins to refine.
		for (int basinId : this.pBasinToRefineIds) {
			subBasinPath = Paths.get(outFolder, String.valueOf(basinId)).toString();
			// run @GeoframeInputBuilder on the basin to zoom.
			pm.message("start work on basin id:" + basinId);
			this.runSmallerBasin(basinId);
			pm.message("create merged file:" + basinId);
			// update the basins of the statistic and the geometries
			// (delete old and insert zoomed)
			this.mergeBasinList(basinId);
			this.mergeFeatures(basinId);
		}
		// copy and move files to new folder
		if (outFolderWithSmallBasin == null || outFolderWithSmallBasin.isEmpty()) {
			outFolderWithSmallBasin = Paths.get(outFolder, "small_basin").toString();
		}
		File newFolderWithSmallBasin = new File(outFolderWithSmallBasin);
		if (!newFolderWithSmallBasin.exists()) {
			newFolderWithSmallBasin.mkdirs();
		}

		dumpVector(newNetworkFC, Paths.get(outFolderWithSmallBasin, NETWORK_SHP_FILE_PATH).toString());
		dumpVector(newBasinFC, Paths.get(outFolderWithSmallBasin, SUBBASINS_COMPLETE_SHP_FILE_PATH).toString());
		this.createBasinCsv();
		// merge topology files.
		this.createTopology();
		this.copyBasin();

	}

	/**
	 * Execute {@link GeoframeInputsBuilder} on the entire basin.
	 * 
	 * @throws Exception
	 */
	private void prepareBigBasin() throws Exception {
		GeoframeInputsBuilder g = new GeoframeInputsBuilder();
		g.inPitfiller = inPitfiller;
		g.inDrain = inDrain;
		g.inTca = inTca;
		g.inNet = inNet;
		g.inSkyview = inSkyview;
		g.inBasins = inBasins;
		g.inLakes = inLakes;
		g.inGeoframeTopology = inGeoframeTopology;
		g.outFolder = outFolder;
		g.process();

	}

	/**
	 * Run MyGeoframeBuilder to refine subbasin.
	 * 
	 * 
	 * @throws Exception
	 */
	private void runSmallerBasin(int basinId) throws Exception {
		HashMap<String, String> inputData = prepareInputData(basinId);
		GeoframeInputsBuilder myGeoframeBuilder = new GeoframeInputsBuilder();
		myGeoframeBuilder.inPitfiller = inputData.get(inPitfiller);
		myGeoframeBuilder.inDrain = inputData.get(inDrain);
		myGeoframeBuilder.inTca = inputData.get(inTca);
		myGeoframeBuilder.inNet = inputData.get(inNet);
		myGeoframeBuilder.inSkyview = inputData.get(inSkyview);
		myGeoframeBuilder.inBasins = inputData.get(inBasins);
		myGeoframeBuilder.inLakes = inputData.get(inLakes);
		myGeoframeBuilder.inGeoframeTopology = inputData.get(inGeoframeTopology);
		myGeoframeBuilder.outFolder = inputData.get(outFolder);
		myGeoframeBuilder.process();

	}

	/*
	 * Copy and move basins folders to clean directory.
	 */
	private void copyBasin() throws IOException {
		String sep = BASIN_CSV_SEPARATOR;
		for (int i = 1; i < bigBasinList.size(); i++) {
			String line = bigBasinList.get(i);
			int id = Integer.valueOf(line.split(sep)[0]);
			File folderToPut = new File(Paths.get(outFolder, String.valueOf(id)).toString());
			if (folderToPut.exists()) {
				String pathId = String.valueOf(id);
				File in = Paths.get(outFolder, pathId).toFile();
				File out = Paths.get(outFolderWithSmallBasin, pathId).toFile();
				FileUtils.copyDirectory(in, out);
			} else {
				for (int j : pBasinToRefineIds) {
					folderToPut = new File(Paths.get(outFolder, String.valueOf(j), String.valueOf(id)).toString());
					if (folderToPut.exists()) {
						String pathBasin = String.valueOf(j);
						String pathId = String.valueOf(id);
						File in = Paths.get(outFolder, pathBasin, pathId).toFile();
						File out = Paths.get(outFolderWithSmallBasin).toFile();
				
						FileUtils.moveDirectoryToDirectory(in, out, true);
					}
				}
			}
		}
	}

	private void createTopology() throws IOException {
		File topoFile = new File(inGeoframeTopology);
		File newTopoFile = new File(topoFile.getParentFile(),
				FileUtilities.getNameWithoutExtention(topoFile) + "_lakes.txt");
		if (newTopoFile.exists()) {
			topoFile = newTopoFile;
		}
		List<String> bigBasinData = FileUtilities.readFileToLinesList(topoFile.toString());
		for (int basinId : pBasinToRefineIds) {
			String subBasinTopo = Paths.get(outFolder, String.valueOf(basinId), TOPOLOGY_CSV_PATH).toString();
			List<String> subBasinData = FileUtilities.readFileToLinesList(subBasinTopo);
			String start = String.valueOf(basinId) + WHITE_SPACE;
			String end = WHITE_SPACE + basinId;
			List<String> newTopo = new ArrayList<String>();
			for (String s : bigBasinData) {

				if (s.startsWith(start)) {
					String inBasin = s.split(WHITE_SPACE)[1];
					for (String f : subBasinData) {
						if (f.endsWith(TOPOLOGY_OUTLET)) {
							f = f.replace(TOPOLOGY_OUTLET, WHITE_SPACE + inBasin);
						}
						newTopo.add(f);
					}
				} else if (s.endsWith(end)) {
					String outBasin = s.split(WHITE_SPACE)[0];

					String inBasin = this.getInBasin(Integer.valueOf(outBasin), newNetworkFC, newBasinFC);
					newTopo.add(outBasin + WHITE_SPACE + inBasin);

				} else {
					newTopo.add(s);
				}
			}
			bigBasinData = newTopo;
		}

//		for(String b : bigBasinData) {
//			if (b.endsWith("null")) {
//				int index = bigBasinData.indexOf(b);
//				String outBasin = b.split(WHITE_SPACE)[0];
//				String inBasin = this.getInBasin(Integer.valueOf(outBasin), newNetworkFC, newBasinFC);
//				b.replace("null",inBasin );
//				bigBasinData.set(index, b);
//
//			}
//		}

		ListIterator<String> it = bigBasinData.listIterator();

		while (it.hasNext()) {
			String b = it.next();

			if (b.endsWith("null")) {
				String outBasin = b.split(WHITE_SPACE)[0];
				String inBasin = this.getInBasin(Integer.valueOf(outBasin), newNetworkFC, newBasinFC);
				b=b.replace("null", inBasin);
				it.set(b);

			}
		}

		File csvFile = new File(Paths.get(outFolderWithSmallBasin).toString(), TOPOLOGY_CSV_PATH);
		if (!csvFile.exists() || doOverWrite) {
			FileUtilities.writeFile(String.join(NEW_LINE, bigBasinData).toString(), csvFile);
		}

	}

	/*
	 * Prepare the input data used to GeoFramInputBuilder to zoom. Store all value
	 * into an HashMap where the key are the input of whole basin.
	 * 
	 * 1) Check if there are lakes or measurements points into the basins. If so,
	 * cut these files on the basin. 2) Crop raster file 3) run NetNumbering on the
	 * basin
	 * 
	 * N.B. if lakes are on the boundary then stop!!!!!! N.B. all files are saved
	 * into the basin folder, inputs are files path.
	 */
	private HashMap<String, String> prepareInputData(int basinId) throws Exception {
		/*
		 * TODO search into the folder the file with extension
		 */
		String fileType = ".asc";
		HashMap<String, String> inputData = new HashMap<String, String>();
		inputData.put(outFolder, subBasinPath);
		inputData.put(inSkyview, getSubbasinLocalPath(subBasinPath, "sky_", fileType, basinId));
		GridCoverage2D maskCoverage = getRaster(this.getSubbasinLocalPath(subBasinPath, "dtm_", fileType, basinId));

		/*
		 * TODO it can be work with input shapefile (without crop????)
		 */
		String boundary = this.getSubbasinLocalPath(subBasinPath, "subbasins_complete_ID_", ".shp", basinId);
		SimpleFeatureCollection basinsFC = getVector(boundary);

		List<SimpleFeature> basinsFList = FeatureUtilities.featureCollectionToList(basinsFC);
		Envelope basinEnvelope = ((Geometry) basinsFList.get(0).getDefaultGeometry()).getEnvelopeInternal();

		if (inLakes != null) {
			SimpleFeatureCollection lakeFC = getVector(inLakes);
			List<SimpleFeature> lakesFList = FeatureUtilities.featureCollectionToList(lakeFC);
			List<SimpleFeature> myLake = new ArrayList<SimpleFeature>();
			SimpleFeature myBasin = basinsFList.get(0);
			double dx = CoverageUtilities.getRegionParamsFromGridCoverage(maskCoverage).getXres();
			for (SimpleFeature l : lakesFList) {

				if (((Geometry) myBasin.getDefaultGeometry())
						.intersects(((Geometry) l.getDefaultGeometry()).buffer(dx))) {
					throw new ModelsRuntimeException("Lake on boundary not yet supported.", this);
				}
				if (((Geometry) myBasin.getDefaultGeometry()).contains((Geometry) l.getDefaultGeometry())) {
					myLake.add(l);
				}
			}
			String inLakesSub = this.getSubbasinLocalPath(subBasinPath, "lake_", ".shp", basinId);
			dumpVector(DataUtilities.collection(myLake), inLakesSub);
			inputData.put(inLakesSub, inLakesSub);
		}

		List<SimpleFeature> myPoint = new ArrayList<SimpleFeature>();
		if (inPoints != null) {
			SimpleFeatureCollection pointsFC = getVector(inPoints);
			List<SimpleFeature> pointsFList = FeatureUtilities.featureCollectionToList(pointsFC);
			myPoint = new ArrayList<SimpleFeature>();
			SimpleFeature myBasin = basinsFList.get(0);
			for (SimpleFeature l : pointsFList) {
				if (((Geometry) myBasin.getDefaultGeometry()).contains((Geometry) l.getDefaultGeometry())) {
					myPoint.add(l);
				}
			}
			String inPointSub = this.getSubbasinLocalPath(subBasinPath, "point_", ".shp", basinId);
			dumpVector(DataUtilities.collection(myPoint), inPointSub);
			inputData.put(inPointSub, inPointSub);
		}
		// crop tca and flow on the basin to zoom.
		String outTca = this.getSubbasinLocalPath(subBasinPath, "tca_bounded_ID", fileType, basinId);
		GridCoverage2D tca = getRaster(inTca);
		GridCoverage2D flow = getRaster(inDrain);

		CoordinateReferenceSystem crs = tca.getCoordinateReferenceSystem();
		ReferencedEnvelope basinRefEnvelope = new ReferencedEnvelope(basinEnvelope, crs);
		GridCoverage2D tcaClipped = CoverageUtilities.clipCoverage(tca, basinRefEnvelope);
		GridCoverage2D flowClipped = CoverageUtilities.clipCoverage(flow, basinRefEnvelope);

		tcaClipped = CoverageUtilities.coverageValuesMapper(tcaClipped, maskCoverage);
		flowClipped = CoverageUtilities.coverageValuesMapper(flowClipped, maskCoverage);

		File tcaFile = new File(outTca);
		if (!tcaFile.exists() || doOverWrite) {
			dumpRaster(tcaClipped, tcaFile.getAbsolutePath());
		}
		inputData.put(inTca, outTca);
//		File flowFile = new File(inFlow);
//		if (!flowFile.exists() || doOverWrite) {
//			dumpRaster(this.addBorderIntNoValue(flowClipped), flowFile.getAbsolutePath());
//		}

		inputData.put(inPitfiller, getSubbasinLocalPath(subBasinPath, "dtm_", fileType, basinId));
		String netFilePath = getSubbasinLocalPath(subBasinPath, "net_", fileType, basinId);
		GridCoverage2D net = this.addBorderIntNoValue(getRaster(netFilePath));
		inputData.put(inNet, netFilePath);
		this.runNN(subBasinPath, flowClipped, net, tcaClipped, fileType, inputData, DataUtilities.collection(myPoint),
				basinId);
		return inputData;

	}

	/**
	 * 
	 * @param path
	 * @param inFlow
	 * @param inNet
	 * @param inTca
	 * @param fileType
	 * @param map
	 * @param inPointsVector
	 * @param basinId
	 * @throws Exception
	 */
	private void runNN(String path, GridCoverage2D inFlow, GridCoverage2D inNet, GridCoverage2D inTca, String fileType,
			HashMap<String, String> map, SimpleFeatureCollection inPointsVector, int basinId) throws Exception {
		String outDesiredBasins = this.getSubbasinLocalPath(path, "nn_basin_desired_", fileType, basinId);
		OmsNetNumbering omsnetnumbering = new OmsNetNumbering();
		GridCoverage2D mo = runMarkOutlet(path, inFlow, basinId, inDrain, map, fileType);
		omsnetnumbering.inFlow = mo;
		omsnetnumbering.inTca = this.addBorderIntNoValue(inTca);
		omsnetnumbering.inNet = inNet;
		omsnetnumbering.inPoints = inPointsVector;
		omsnetnumbering.pDesiredArea = pDesiredArea;
		omsnetnumbering.pDesiredAreaDelta = pDesiredAreaDelta;
		omsnetnumbering.doProcess = doProcess;
		omsnetnumbering.doReset = doReset;
		omsnetnumbering.process();
		int valueToAdd = Integer.valueOf(String.valueOf(basinId) + GeoFrameInputBuilderWithSubBasin.basinPathSeparator);
		String[] topoBasin = omsnetnumbering.outGeoframeTopology.split(WHITE_SPACE + "|" + NEW_LINE);

		StringBuilder newTopology = new StringBuilder();
		for (int i = 0; i < topoBasin.length; i = i + 2) {
			int firstValue = valueToAdd + Integer.valueOf(topoBasin[i]);

			int secondValueValue = Integer.valueOf(topoBasin[i + 1]);
			if (secondValueValue != 0) {
				secondValueValue = secondValueValue + valueToAdd;
			}
			newTopology.append(firstValue + WHITE_SPACE + secondValueValue + NEW_LINE);

		}
		File csvFile = new File(Paths.get(outFolder, String.valueOf(basinId)).toString(), TOPOLOGY_CSV_PATH);
		if (!csvFile.exists() || doOverWrite) {
			FileUtilities.writeFile(newTopology.toString(), csvFile);
		}
		map.put(this.inGeoframeTopology, csvFile.getAbsolutePath());
		RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(omsnetnumbering.inTca);
		int nCols = regionMap.getCols();
		int nRows = regionMap.getRows();
		WritableRaster basinWR = CoverageUtilities
				.renderedImage2WritableRaster(omsnetnumbering.outDesiredBasins.getRenderedImage(), false);
		WritableRandomIter basinIter = RandomIterFactory.createWritable(basinWR, null);
		for (int r = 0; r < nRows; r++) {
			for (int c = 0; c < nCols; c++) {
				int basinNumber = basinIter.getSample(c, r, 0);
				if (!HMConstants.isNovalue(basinNumber))
					basinWR.setSample(c, r, 0, basinIter.getSample(c, r, 0) + valueToAdd);
			}
		}
		GridCoverage2D desiredBasin = removeBorderIntNoValue(CoverageUtilities.buildCoverage("", basinWR, regionMap,
				omsnetnumbering.inTca.getCoordinateReferenceSystem()));

		dumpRaster(desiredBasin, outDesiredBasins);

		if (omsnetnumbering.outGeoframeTopology != null && omsnetnumbering.outGeoframeTopology.trim().length() > 0) {
			String outGeoframeTopology = outFolder + File.pathSeparator + String.valueOf(basinId) + File.pathSeparator
					+ TOPOLOGY_CSV_PATH;
			FileUtilities.writeFile(omsnetnumbering.outGeoframeTopology, new File(outGeoframeTopology));
			map.put(inGeoframeTopology, omsnetnumbering.outGeoframeTopology);

		}
		map.put(inBasins, outDesiredBasins);

	}

	/*
	 * add no value on the raster edge
	 * 
	 */
	private GridCoverage2D addBorderIntNoValue(GridCoverage2D rasterToAddBorder) throws Exception {
		RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(rasterToAddBorder);
		int nCols = regionMap.getCols() + 2;
		int nRows = regionMap.getRows() + 2;
		WritableRaster boundedWR = CoverageUtilities.createWritableRaster(nCols, nRows, Integer.class, null,
				HMConstants.intNovalue);
		RegionMap newRegionMap = regionMap.toSubRegion(regionMap.getNorth() + regionMap.getYres(),
				regionMap.getSouth() - regionMap.getYres(), regionMap.getWest() - regionMap.getXres(),
				regionMap.getEast() + regionMap.getXres());
		WritableRaster toBoundedWR = CoverageUtilities
				.renderedImage2WritableRaster(rasterToAddBorder.getRenderedImage(), false);
		WritableRandomIter toBoundedIter = RandomIterFactory.createWritable(toBoundedWR, null);
		for (int r = 1; r < nRows - 1; r++) {
			for (int c = 1; c < nCols - 1; c++) {
				boundedWR.setSample(c, r, 0, toBoundedIter.getSample(c - 1, r - 1, 0));
			}
		}
		return CoverageUtilities.buildCoverage("markoutlets", boundedWR, newRegionMap,
				rasterToAddBorder.getCoordinateReferenceSystem());
	}

	/*
	 * remove no value edge
	 * 
	 */
	private GridCoverage2D removeBorderIntNoValue(GridCoverage2D rasterToRemoveBorder) throws Exception {
		RegionMap regionMap = CoverageUtilities.getRegionParamsFromGridCoverage(rasterToRemoveBorder);
		int nCols = regionMap.getCols() - 2;
		int nRows = regionMap.getRows() - 2;
		WritableRaster boundedWR = CoverageUtilities.createWritableRaster(nCols, nRows, Integer.class, null,
				HMConstants.intNovalue);
		RegionMap newRegionMap = regionMap.toSubRegion(regionMap.getNorth() - regionMap.getYres(),
				regionMap.getSouth() + regionMap.getYres(), regionMap.getWest() + regionMap.getXres(),
				regionMap.getEast() - regionMap.getXres());
		WritableRaster toBoundedWR = CoverageUtilities
				.renderedImage2WritableRaster(rasterToRemoveBorder.getRenderedImage(), false);
		WritableRandomIter toBoundedIter = RandomIterFactory.createWritable(toBoundedWR, null);
		for (int r = 0; r < nRows; r++) {
			for (int c = 0; c < nCols; c++) {
				boundedWR.setSample(c, r, 0, toBoundedIter.getSample(c + 1, r + 1, 0));
			}
		}
		return CoverageUtilities.buildCoverage("markoutlets", boundedWR, newRegionMap,
				rasterToRemoveBorder.getCoordinateReferenceSystem());

	}

	private GridCoverage2D runMarkOutlet(String path, GridCoverage2D inFlow, int basinId, String key,
			HashMap<String, String> map, String fileType) throws Exception {
		outMarkOutlet = this.getSubbasinLocalPath(path, "mo_", fileType, basinId);
		OmsMarkoutlets omsmarkoutlets = new OmsMarkoutlets();
		omsmarkoutlets.inFlow = inFlow;
		omsmarkoutlets.pm = pm;
		omsmarkoutlets.doProcess = doProcess;
		omsmarkoutlets.doReset = doReset;
		omsmarkoutlets.process();
		dumpRaster(omsmarkoutlets.outFlow, outMarkOutlet);
		map.put(key, outMarkOutlet);
		// add a novalue border to he result, otherwise netnumbering not works.
		return addBorderIntNoValue(omsmarkoutlets.outFlow);
	}

	/*
	 * get the file path in basin subdirectory
	 */
	private String getSubbasinLocalPath(String path, String prefix, String fileType, int basinId) {
		StringBuilder sb = new StringBuilder(path);
		sb.append(File.separator).append(prefix).append(basinId).append(fileType);
		return sb.toString();
	}

	/**
	 * Merge the subbasin.csv files
	 * 
	 * @throws Exception
	 * @throws ClassNotFoundException
	 */
	private void mergeFeatures(int basinId) throws ClassNotFoundException, Exception {
		String netSubBasin = Paths.get(outFolder, String.valueOf(basinId), NETWORK_SHP_FILE_PATH).toString();
		String boundary = this.getSubbasinLocalPath(Paths.get(outFolder, String.valueOf(basinId)).toString(),
				"subbasins_complete_ID_", ".shp", basinId);
		SimpleFeatureCollection basinsFC = getVector(boundary);

		List<SimpleFeature> basinsFList = FeatureUtilities.featureCollectionToList(basinsFC);
		SimpleFeature myBasin = basinsFList.get(0);
		SimpleFeatureCollection netSubFC = getVector(netSubBasin);
		SimpleFeatureIterator it = newNetworkFC.features();

		List<SimpleFeature> netFList = FeatureUtilities.featureCollectionToList(newNetworkFC);
		List<SimpleFeature> netSubFList = FeatureUtilities.featureCollectionToList(netSubFC);

		while (it.hasNext()) {
			SimpleFeature l = it.next();
			if ((int) l.getAttribute(BASINID_KEY) == basinId) {
				if (((Geometry) myBasin.getDefaultGeometry()).intersects(((Geometry) l.getDefaultGeometry()))) {
					SimpleFeatureIterator itSub = netSubFC.features();
					while (itSub.hasNext()) {
						SimpleFeature inNet = itSub.next();
						if (((Geometry) inNet.getDefaultGeometry()).within(((Geometry) l.getDefaultGeometry()))) {
							netSubFList.remove(inNet);
							inNet.setDefaultGeometry(
									((Geometry) inNet.getDefaultGeometry()).union(((Geometry) l.getDefaultGeometry())));
							netSubFList.add(inNet);
						}
					}

					netFList.remove(l);
				}
			}
		}
		UnionFeatureCollection unionFeatureCollection = new UnionFeatureCollection();
		newNetworkFC = unionFeatureCollection.execute(DataUtilities.collection(netFList),
				DataUtilities.collection(netSubFList));

		String idSubBsin = Paths.get(outFolder, String.valueOf(basinId), SUBBASINS_COMPLETE_SHP_FILE_PATH).toString();
		SimpleFeatureCollection refBasinsFC = getVector(idSubBsin);

		it = newBasinFC.features();

		List<SimpleFeature> basinFList = FeatureUtilities.featureCollectionToList(newBasinFC);
		while (it.hasNext()) {
			SimpleFeature l = it.next();
			if ((int) l.getAttribute(BASINID_KEY) == basinId) {
				basinFList.remove(l);
			}
		}
		newBasinFC = unionFeatureCollection.execute(DataUtilities.collection(basinFList), refBasinsFC);

	}

	/*
	 * To reconstruct topology, get the id of the basin that "flows" into the basin
	 * with id==basinId
	 */
	private String getInBasin(Integer basinId, SimpleFeatureCollection netFC, SimpleFeatureCollection refBasinsFC) {
		// TODO Auto-generated method stub
		SimpleFeatureIterator it = netFC.features();

		while (it.hasNext()) {
			SimpleFeature l = it.next();
			if ((int) l.getAttribute(BASINID_KEY) == basinId) {
				SimpleFeatureIterator itSub = refBasinsFC.features();
				while (itSub.hasNext()) {
					SimpleFeature inBasin = itSub.next();
					if ((int) inBasin.getAttribute(BASINID_KEY) != basinId && ((Geometry) inBasin.getDefaultGeometry())
							.intersects(((Geometry) l.getDefaultGeometry()))) {
						return String.valueOf(inBasin.getAttribute(BASINID_KEY));
					}
				}
			}
		}

		return null;
	}

	/**
	 * ˇ Merge the sub-basin statistic ( remove the basin to zoom and add zoomed
	 * sub-basins)
	 * 
	 * @throws IOException
	 */
	private void mergeBasinList(int basinId) throws IOException {
		String csvID = Paths.get(outFolder, String.valueOf(basinId), SUBBASINS_CSV_PATH).toString();
		String sep = BASIN_CSV_SEPARATOR;
		for (int i = 1; i < bigBasinList.size(); i++) {
			String line = bigBasinList.get(i);
			if (Integer.valueOf(line.split(sep)[0]) == basinId) {
				List<String> newText = FileUtilities.readFileToLinesList(csvID);
				newText = newText.subList(1, newText.size());
				bigBasinList.remove(line);
				bigBasinList.addAll(i, newText);
				return;
			}
		}

	}

	/*
	 * create the basin sat csv file.
	 */
	private void createBasinCsv() throws IOException {
		Iterator<String> iterator = bigBasinList.iterator();
		// read header
		StringBuilder csvText = new StringBuilder();

		if (iterator.hasNext()) {
			String line = iterator.next();
			csvText.append(line);
			csvText.append(NEW_LINE);
		}

		while (iterator.hasNext()) {
			String line = iterator.next();
			csvText.append(line);
			csvText.append(NEW_LINE);

		}
		File csvFile = new File(Paths.get(outFolderWithSmallBasin).toString(), SUBBASINS_CSV_PATH);
		if (!csvFile.exists() || doOverWrite) {
			FileUtilities.writeFile(csvText.toString(), csvFile);
		}

	}

	public static void main(String[] args) throws Exception {
		String path = "/home/andreisd/Documents/project/GWS2022/OMS_project_val_di_non/OMS_project/data/Trentino/Noce/geomorphology/Noce/net100/";
		String pit = path + "Noce_pitfiller.tif";
		String drain = path + "Noce_drain_dir.tif";
		String tca = path + "Noce_Tca.tif";
		String net = path + "Noce_network_100.tif";
		String sky = path + "Noce_skyview.tif";
		String basins = path + "Noce_subbasins_desired.tif";
		String topology = path + "topology_10km.csv";
		String pathData = "/home/andreisd/Documents/project/GWS2022/OMS_project_val_di_non/OMS_project/data/Trentino/Noce/";
		String lakes = pathData + "laghi.shp";
		String point = pathData + "idrometri_netnumbering_network_100_sgiustinacorto.shp";
		String outfolder = path + "km1/geoframe";
		GeoFrameInputBuilderWithSubBasin g = new GeoFrameInputBuilderWithSubBasin();
		g.inPitfiller = pit;
		g.inDrain = drain;
		g.inTca = tca;
		g.inNet = net;
		g.inSkyview = sky;
		g.inBasins = basins;
		g.inLakes = lakes;
		g.inGeoframeTopology = topology;
		g.outFolder = outfolder;
		g.pBasinToRefineIds = new int[] { 4645, 3505, 4506, 4829, 5208 };
		g.pDesiredArea = 1000000.0;
		g.pDesiredAreaDelta = 20.0;
		g.inPoints = point;
		g.pGeoFrameAlreadyDone = true;
		g.process();
	}
}
