//////////////////////////////////////////// 
// See Murray et al. (Nature, 2019) for further information.
// Developed in Google Earth Engine:
// https://code.earthengine.google.com
// Cite:
// Murray N. J., Phinn S. R., DeWitt M., Ferrari R., Johnston R., Lyons M. B., Clinton N., Thau D. & 
// Fuller R. A. (2019) The global distribution and trajectory of tidal flats. *Nature*. 565, 222-225.
//////////////////////////////////////////// 

// 0. Global Variables
var site = ee.Geometry.Polygon([-180, 60, 0, 60, 180, 60, 180, -60, 10, -60, -180, -60], null, false);
var globOptions = { 
  versionID: '_SR',
  outFolder: 'SR',
  startDate: '2014-01-01',
  endDate: '2016-12-31',
  bandSelect: ['green', 'swir1', 'swir2', 'nir', 'red'],
  bands8: ['B3', 'B6', 'B7', 'B5', 'B4'],
  bands7: ['B2', 'B5', 'B7', 'B4', 'B3'], 
  maskAltitude: 100,  
  maskDepth: -100, 
  maskDistance: 50000,
  maskApplySRTM: false,
  parallelScale: 8,
  trainingValidationRatio: 0.0001 
  nTrees: 10, 
  outScale: 30, 
  conPixels: 100
};

// 1. Functions
var landsatFunctions = {
  applyFMask: function(image) {
    // Mask out SHADOW, SNOW, and CLOUD classes. SR data.
    return image
      .updateMask(image.select('cfmask')
      .lt(2)); 
  },

  applyNDWI: function(image) {
    // apply NDWI to an image
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']);
  },

  applyMNDWI: function(image) {
    // apply MNDWI to an image
    var mndwi = image.normalizedDifference(["green","swir1"]);
    return image.select([0], ['mndwi']);
  },

  applyAWEI: function(image) {
    // apply AWEI to an image
    var awei = image.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))");
    return awei.select([0], ['awei']);
  },

  applyNDVI: function(image) {
    // apply NDVI to an image
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']);
  },
};

var reducer = ee.Reducer.min()
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.median(), '', true)
    //.combine(ee.Reducer.count(), '', true) 
    .combine(ee.Reducer.percentile([10, 25, 50, 75,90]), '', true)
    .combine(ee.Reducer.intervalMean(0, 10).setOutputs(['intMn0010']), '', true)
    .combine(ee.Reducer.intervalMean(10, 25).setOutputs(['intMn1025']), '', true)
    .combine(ee.Reducer.intervalMean(25, 50).setOutputs(['intMn2550']), '', true)
    .combine(ee.Reducer.intervalMean(50, 75).setOutputs(['intMn5075']), '', true)
    .combine(ee.Reducer.intervalMean(75, 90).setOutputs(['intMn7590']), '', true)
    .combine(ee.Reducer.intervalMean(90, 100).setOutputs(['intMn90100']), '', true)
    .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090']), '', true)
    .combine(ee.Reducer.intervalMean(25, 75).setOutputs(['intMn2575']), '', true);

// 2. Data Imports & Processing
// vectors
var globCoast = ee.FeatureCollection('ft:1Hsoe_WwULJ23Nuj1wikGzfH_WQMtpDWOR3XpWkHk');
var randomPointsPreComputed = ee.FeatureCollection('ft:1hVC5uIlWZQxtsapNsr5AzcLm1Vzo4I_DqjfNgNmN'); // Precomputed Training

// images
function generateLandsatCollection(){
  var L4collection = ee.ImageCollection('LANDSAT/LT4_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L5collection = ee.ImageCollection('LANDSAT/LT5_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L7collection = ee.ImageCollection('LANDSAT/LE7_SR')
      .filterDate(globOptions.startDate,globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L8collection = ee.ImageCollection('LANDSAT/LC8_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands8, globOptions.bandSelect);
  var collectionFull = ee.ImageCollection(L4collection
      .merge(L5collection)
      .merge(L7collection)
      .merge(L8collection))
      //.filterBounds(site)
      .filter(ee.Filter.intersects('.geo', globCoast.geometry(), null, null, 1000))
      .filterMetadata('WRS_ROW', 'less_than', 120); 
  return collectionFull;
}
var collection = generateLandsatCollection();

// Data processing
var covariates = {
    aweiReduced: collection.map(landsatFunctions.applyAWEI)
        .reduce(reducer, globOptions.parallelScale),
    ndwiReduced: collection.map(landsatFunctions.applyNDWI)
        .reduce(reducer, globOptions.parallelScale),
    mndwiReduced: collection.map(landsatFunctions.applyMNDWI)
        .reduce(reducer, globOptions.parallelScale),
    ndvi: collection.map(landsatFunctions.applyNDVI)
        .reduce(ee.Reducer.intervalMean(10, 90)
        .setOutputs(['intMn1090'])),
    nirBand: collection.select(['nir'])
        .reduce(ee.Reducer.intervalMean(10, 90)
        .setOutputs(['intMn1090'])),
    swir1Band: collection.select(['swir1'])
        .reduce(ee.Reducer.intervalMean(10, 90)
        .setOutputs(['intMn1090'])),
    etopo: ee.Image('NOAA/NGDC/ETOPO1')
        .select(['bedrock'], ['etopo'])
        .resample('bicubic'),
    swOccurrence: ee.Image('JRC/GSW1_0/GlobalSurfaceWater')
        .select(['occurrence'], ['surfaceWater'])
        .unmask()
};
var trainComposite = covariates.aweiReduced
    .addBands(covariates.ndwiReduced)
    .addBands(covariates.mndwiReduced)
    .addBands(covariates.ndvi)
    .addBands(covariates.nirBand)
    .addBands(covariates.swir1Band)
    .addBands(covariates.etopo)
    .addBands(covariates.swOccurrence);
var bands = trainComposite.bandNames();

// 3. Masking
var coastMask = globCoast
    .distance(globOptions.maskDistance).gte(-20); 
var topoMask = covariates.etopo
    .gte(globOptions.maskDepth)
    .and(covariates.etopo.lte(globOptions.maskAltitude));
var topoMask = topoMask.updateMask(topoMask);
if (globOptions.maskApplySRTM) {
  var srtm = ee.Image('USGS/SRTMGL1_003')
    .lte(0);
  var finalMask = coastMask.multiply(topoMask)
    .rename('datamask')
    .byte()
    .updateMask(srtm);
} else {
  var finalMask = coastMask.multiply(topoMask)
    .rename('datamask')
    .byte();
}

// 4. Training Data
var predictorSet = randomPointsPreComputed;
var predictorSubSet = predictorSet
  .filter(ee.Filter.neq('ndwi_stdDev', null)) 
  .randomColumn('random', 0);
var trainingSet = predictorSubSet
  .filter(ee.Filter.gte('random',globOptions.trainingValidationRatio));
var validationSet = predictorSubSet
  .filter(ee.Filter.lt('random',globOptions.trainingValidationRatio));

// 5. Classify 
var classifier = ee.Classifier.randomForest({
    numberOfTrees: globOptions.nTrees, 
    variablesPerSplit: 0,
    bagFraction: 0.5,
    seed: 0})
  .train(trainingSet, 'CLASS', bands)
  .setOutputMode('CLASSIFICATION');
var classified = trainComposite.select(bands)
  .classify(classifier);
var finalOut = classified.byte()
  .mask(finalMask); 

// Postprocess
 var finalOut = finalOut.mask(finalOut
  .connectedPixelCount(globOptions.conPixels)
  .gte(globOptions.conPixels)); 
 
var finalOut = finalOut.updateMask(finalOut.eq(2)); // tf only


// Extra post-process
var terrestrial5k = ee.FeatureCollection("users/murrnick/tidalFlat/postProcessing/simBoundary_Buff5k");
var coastMask5k = ee.Image(1).clip(terrestrial5k);
var invcoastMask5k = coastMask5k.mask(coastMask5k.mask().not());
var finalOut = finalOut.updateMask(invcoastMask5k);

// 8. Run
print(ee.Serializer.toReadableJSON(finalOut));