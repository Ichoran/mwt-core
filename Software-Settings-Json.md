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

The top-level struture is a JSON object with up to five subcategories:

1. `"segmentation"` - This key is associated with a JSON object that lists
the image segmentation parameters that form the bulk of the settings in MWT-core.

2. `"output"` - This key is associated with a JSON object that says where and
what is being produced, including how long the recording was supposed to be.

3. `"masks"` - This key is associated with a JSON object that specifies how the
regions of analysis should be restricted.

4. `"stimuli"` - Although MWT-core does not itself deliver stimuli, it needs
to know when they occur to mark the output file.  This key is associated with
an array of JSON objects that specify what stimuli are expected.

5. `"custom"` - Any front-end specific values go under this key.

In addition, the `"software"` key may be associated with an identifying string that
specifies the software that took the data, e.g. `"MWT LabView 1.5.0-M2"`

## The `segmentation` object

The segmentation object contains the following key-value pairs to describe
how worms are detected and separated from background.

1. `"dark"` : `true` or `false`.  Specifies whether objects are brighter or darker than background.
2. `"binning"`: Integer number specifying how much to bin.  Assume a default of 1.
3. `"contrast"`: number indicating percentage of full dynamic range darker/brighter than background an object should be.
4. `"contrast-hysteresis"`: number indicating fraction of contrast to fill to (closer to background)
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
6. `"snapshots"`: a number indicating how frequently images are saved.
7. `"dbde"`: `true` or `false` depending on whether images are saved to a DBDE file.  (If false, it's separate PNGs.)

## The `masks` object

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

## The `"custom"` value

MWT-LabView uses two fields in an object here.

1. `"bit-depth"` - Either a number, which is the bit depth, or the string `"camera"` meaning that the bit depth is read from the camera
2. `"auto-start"` - A number indicating the time in seconds after which to automatically begin recording
