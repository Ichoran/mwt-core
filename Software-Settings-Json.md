# Saving Multi-Worm Tracker settings to JSON

It's useful to have a standard form for saving tracker settings.  JSON is
a widely-used format for settings because of its simplicity and because
there are many readers and writers for JSON.

To read more about the JSON format, see [json.org](http://json.org).

This document recommends a format for saving Multi-Worm Tracker settings.

The MWT-LabView front-end uses this format.

## Overall structure

In general, all entries are optional.  If a key-value pair is missing, some
sensible default should be used.

All numbers should be finite.

The top-level struture is a JSON object with up to eight subcategories:

1. `"timestamp"` - This key is associated with a string containing an ISO-8601 timestamp.  Only local and UTC timezones are supported.

1. `"software"` - This key is associated with a string value that specifies the version of the software that created the file.

1. `"segmentation"` - This key is associated with a JSON object that lists
the image segmentation parameters that form the bulk of the settings in MWT-core.

2. `"output"` - This key is associated with a JSON object that says where and
what is being produced, including how long the recording was supposed to be.

3. `"masks"` - This key is associated with a JSON array that specifies how the
regions of analysis should be restricted.

4. `"stimuli"` - Although MWT-core does not itself deliver stimuli, it needs
to know when they occur to mark the output file.  This key is associated with
an array of JSON objects that specify what stimuli are expected.

5. `"reference"` - Locations of reference objects.

5. `"custom"` - Any front-end specific values go under this key.


## The `segmentation` object

The segmentation object contains the following key-value pairs to describe
how worms are detected and separated from background.

1. `"dark"` : `true` or `false`.  Specifies whether objects are brighter or darker than background.
2. `"binning"`: Integer number specifying how much to bin.  Assume a default of 1.
3. `"contrast"`: number indicating fraction of full dynamic range darker/brighter than background an object should be.
4. `"contrast-hysteresis"`: number indicating fraction of contrast less than limit to fill to (closer to background)
5. `"size-min"`: number of pixels indicating the smallest allowable object
6. `"size-max"`: number of pixels indicating the largest allowable object
7. `"size-hysteresis"`: fraction above/below max/min you can go with an object you're already tracking before you discard it
8. `"alpha"`: a number indicating negative base two exponent of the adaptation rate (implemented as a bit shift) of the background.
9. `"bands"`: a number indicating how many bands to subdivide the image into for ongoing background adaptation and object detection.
10. `"border"`: a number indicating how many pixels around a detected object are scanned in the next frame.
11. `"divisive"`: `true` or `false` depending on whether divisive normalization is used (default is false, meaning subtractive).

## The `output` object

1. `"prefix"`: a string containing the filename prefix / experiment description.
2. `"tracker"`: a string containing an identifier for the specific tracker used.
3. `"outline"`: `true` or `false` depending on whether outlines are generated.  I don't know why you'd ever want this false!
4. `"skeleton"`: `true` or `false` depending on whether skeletons are generated.
5. `"duration"`: a number indicating the expected duration of the experiment in seconds.
6. `"snapshots"`: a number indicating how frequently images are saved in seconds between images (0 indicates no snapshots)
7. `"dbde"`: `true` or `false` depending on whether images are saved to a DBDE file.  (If false, it's separate PNGs.)

## The `masks` array

The `masks` object contains an array of shape-action objects.  The first one is assumed to be set while the rest are adds or cuts.

The shape-action objects contain the following:

1. `"ellipse"`: `true` or `false`; if false `"shape"` is a rectangular mask object, otherwise it's elliptical (see below)
2. `"shape"`: an array of four integers describing either a rectangular mask or an elliptical mask depending on the value of `"ellipse"`
3. `"cut"`: `true` or `false`; if false, the new shape is added to the existing mask.  Otherwise, it is cut.  The default should be false.

### Rectangular mask values

A rectangular mask is specified by its two corners.  The four elements of the array are interpreted as:

1. `x1`, the integer value of the smaller x coordinate
2. `y1`, the integer value of the smaller y coordinate
3. `x2`, the integer value of the larger x coordinate
4. `y2`, the integer value of the larger y coordinate

### Elliptical mask values

An elliptical mask is specified by a center and radii.  The four elements of the array are interpreted as:

1. `cx`, the integer value of the x coordinate of the center point of the ellipse
2. `cy`, the integer value of the y coordinate of the center point of the ellipse
3. `rx`, the integer value of the x radius of the ellipse
4. `ry`, the integer value of the y radius of the ellipse

## The `stimuli` array

Stimuli are specified in an array.  Each entry is a JSON object with the following key-value pairs:

1. `"device"` - A string describing the type of stimulus device.  Usually this will be `"Ticklish"`
2. `"id"` - A string that specifies the particular identity of the device.
3. `"port"` - A string that describes a port or other access information for the device (as of when the program was run)
4. `"channel"` - A string that describes the output channel (for Ticklish, it's a single letter from `"A"` to `"X"`).
5. `"description"` - A string that describes the type of stimulus (corresponding to old-style tap/puff etc.)
6. `"delay"` - A non-negative number that specifies the time in seconds before the first stimulus of this type
7. `"interval"` - A non-negative number that specifies the time between stimuli
8. `"count"` - A non-negative number that specifies the number of stimuli
9. `"high"` - How long the stimulus is active
10. `"pulses"` - A number that, if greater than 1, indicates how many repeats of high stimulus are given in a row
11. `"pulse-interval"` - Time between pulses within one stimulus (zero is fine if there is only a single high)
12. `"shape"` - One or two letters expressing the shape; for now, probably always `"u"` for upright digital pulses
13. `"more"` - A nested stimulus object that specifies what happens immediately after the first stimulus train is done
14. `"custom"` - A JSON value that contains any custom settings to the stimulator.  For Ticklish, this is the command-string sent to the device to set up the stimulus.

## The `"reference"` object

Reference objects are specified in an object that consists primarily of an array of x,y coordinates.  The object contains three fields:

1. `"low-intensity"`: this key is associated with a number specifying the lower bound of the threshold for a reference object

2. `"high-intensity"`: this key is associated with a number specifying the upper bound of the threshold for a reference objects

3. `"coordinates"`: this is an array of objects; each object has keys `"x"` and `"y"` whose values are numbers that specify the coordinates of the center of the object.

## The `"custom"` value

MWT-LabView uses two fields in an object here.

1. `"bit-depth"` - Either a positive integer, which is the bit depth, or 0, meaning that the bit depth is read from the camera
2. `"auto-start"` - A number indicating the time in seconds after which to automatically begin recording

## An example

```
{
  "timestamp": "2016-12-19T14:42:30",
  "software": "Rex Kerr hand-generated",
  "segmentation": {
    "dark": true,
    "binning": 1,
    "contrast": 0.1,
    "contrast-hysteresis": 0.4,
    "size-min": 50,
    "size-max": 1000,
    "size-hysteresis": 0.1,
    "alpha": 5,
    "bands": 10,
    "border": 20,
    "divisive": false
  },
  "output": {
    "prefix": "N2_72h_01a",
    "tracker": "green",
    "outline": true,
    "skeleton": false,
    "duration": 125.5,
    "snapshots": 0,
    "dbde": false
  },
  "masks": [
    { "ellipse": false, "shape": [100, 100, 1948, 1948] }
  ],
  "stimuli": [
    {
      "device": "Ticklish",
      "id": "mgbv",
      "port": "COM3",
      "channel": "A",
      "description": "tap",
      "delay": 60,
      "interval": 5,
      "count": 13,
      "high": 0.05,
      "pulses": 1,
      "shape": "u",
      "custom": "A=120.0500;60.0000;0.050000;4.950000;0.050000;0.000001u"
    }
  ],
  "reference": {
    "low-intensity": 0,
    "high-intensity": 16000,
    "coordinates": [
      { "x": 1235, "y": 1441 },
      { "x": 551, "y": 1153 }
    ]
  },
  "custom": {
    "bit-depth": 0,
    "auto-start": 2
  }
}
```
