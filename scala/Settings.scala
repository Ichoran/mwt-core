package lab.kerrr.mwt.settings

import kse.jsonal._
import kse.jsonal.JsonConverters._
import kse.flow._
import kse.maths._

object ImplicitMwtSettingsJson {
  implicit val implicitMwtSettingsSegmentation: FromJson[Segmentation] = Segmentation
  implicit val implicitMwtSettingsOutput: FromJson[Output] = Output
  implicit val implicitMwtSettingsMask: FromJson[Mask] = Mask
  implicit val implicitMwtSettingsStimulus: FromJson[Stimulus] = Stimulus
  implicit val implicitMwtSettingsCoordinate: FromJson[Coordinate] = Coordinate
  implicit val implicitMwtSettingsReferences: FromJson[References] = References
  implicit val implicitMwtSettingsCustomLabView: FromJson[CustomLabView] = CustomLabView
  implicit val implicitMwtSettingsSettings: FromJson[Settings] = Settings
}

case class Segmentation(
  dark: Boolean, binning: Int,
  contrast: Double, contrastHyst: Double,
  sizeMin: Int, sizeMax: Int, sizeHyst: Double,
  alpha: Int, bands: Int, border: Int,
  divisive: Boolean
) extends AsJson {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def json = Json ~ 
    ("dark", dark) ~ ("binning", binning max 1) ~
    ("contrast", sigfig(contrast)) ~ ("contrast-hysteresis", sigfig(contrastHyst)) ~
    ("size-min", sizeMin max 1) ~ ("size-max", sizeMax max sizeMin max 1) ~ ("size-hysteresis", sigfig(sizeHyst)) ~
    ("alpha", alpha.clip(1, 15)) ~ ("bands", bands max 1) ~ ("border", border max 0) ~
    ("divisive", divisive) ~
    Json
}
object Segmentation extends FromJson[Segmentation] {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  val default = new Segmentation(true, 1, 0.1f, 0.4f, 50, 2000, 0.1f, 5, 10, 20, false)
  def parse(j: Json): Either[JastError, Segmentation] = j match {
    case o: Json.Obj =>
      val dark = if (o contains "dark")                o("dark").to[Boolean].OUT                          else default.dark
      val bin  = if (o contains "binning")             o("binning").to[Int].OUT                           else default.binning
      val cont = if (o contains "contrast")            o("contrast").to[Double].OUT.fn(sigfig)            else default.contrast
      val cHys = if (o contains "contrast-hysteresis") o("contrast-hysteresis").to[Double].OUT.fn(sigfig) else default.contrastHyst
      val sMin = if (o contains "size-min")            o("size-min").to[Int].OUT                          else default.sizeMin
      val sMax = if (o contains "size-max")            o("size-max").to[Int].OUT                          else default.sizeMax
      val sHys = if (o contains "size-hysteresis")     o("size-hysteresis").to[Double].OUT.fn(sigfig)     else default.sizeHyst
      val alph = if (o contains "alpha")               o("alpha").to[Int].OUT                             else default.alpha
      val band = if (o contains "bands")               o("bands").to[Int].OUT                             else default.bands
      val bord = if (o contains "border")              o("border").to[Int].OUT                            else default.border
      val divs = if (o contains "divisive")            o("divisive").to[Boolean].OUT                      else default.divisive
      if (bin < 0 || bin > 9999) return Left(JastError(f"Absurd binning value: $bin"))
      if (cont <= 0 || cont >= 1.0 || !cont.finite) return Left(JastError(f"Absurd contrast fraction: $cont"))
      if (cHys < 0 || !cHys.finite) return Left(JastError(f"Absurd contrast hysteresis value: $cHys"))
      if (sMin < 0) return Left(JastError(f"Absurd minimim size: $sMin"))
      if (sMax < sMin) return Left(JastError(f"Absurd maximum size: $sMax (min is $sMin)"))
      if (alph < 1 || alph > 15) return Left(JastError(f"Invalid alpha value; 1-15 is okay, but got $alph"))
      if (band < 0) return Left(JastError(f"Absurd number of bands for detection of new worms: $band"))
      if (bord < 0) return Left(JastError(f"Absurd number of border pixels: $bord"))
      Right(new Segmentation(dark, bin max 1, cont, cHys, sMin, sMax, sHys, alph, band, bord, divs))
    case Json.Null   => Right(default)
    case _           => Left(JastError("Expected JSON object for Segmentation Parameters but did not get one"))
  }
}

case class Output(
  prefix: String, tracker: String,
  outline: Boolean, skeleton: Boolean,
  duration: Double, snapshots: Double, dbde: Boolean
) extends AsJson {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def json = Json ~
    ("prefix", if (prefix.isEmpty) "_" else prefix) ~ ("tracker", tracker) ~
    ("outline", outline) ~ ("skeleton", skeleton) ~
    ("duration", sigfig(duration)) ~ ("snapshots", sigfig(snapshots)) ~ ("dbde", dbde) ~
    Json
}
object Output extends FromJson[Output]{
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  val default = new Output("", "", true, true, 1, 0, false)
  def parse(j: Json): Either[JastError, Output] = j match {
    case o: Json.Obj =>
      val pfix = if (o contains "prefix")    o("prefix").to[String].OUT               else default.prefix
      val targ = if (o contains "tracker")   o("tracker").to[String].OUT              else default.tracker
      val outl = if (o contains "outline")   o("outline").to[Boolean].OUT             else default.outline
      val skel = if (o contains "skeleton")  o("skeleton").to[Boolean].OUT            else default.skeleton
      val dura = if (o contains "duration")  o("duration").to[Double].OUT.fn(sigfig)  else default.duration
      val snap = if (o contains "snapshots") o("snapshots").to[Double].OUT.fn(sigfig) else default.snapshots
      val dbde = if (o contains "dbde")      o("dbde").to[Boolean].OUT                else default.dbde
      if (dura < 0 || !dura.finite) return Left(JastError(f"Absurd duration, $dura"))
      if (snap < 0 || !snap.finite) return Left(JastError(f"Absurd snapshots value, $snap"))
      Right(new Output(pfix, targ, outl, skel, dura, snap, dbde))
    case Json.Null   => Right(default)
    case _           => Left(JastError("Expected JSON object for Output but did not get one."))
  }
}

case class Mask(ellipse: Boolean, shape: Array[Int], cut: Boolean) extends AsJson {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def json = Json ~ ("ellipse", ellipse) ~ ("shape", Json(shape)) ~ ("cut", cut) ~ Json 
}
object Mask extends FromJson[Mask] {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def default = new Mask(false, Array(0, 0, 0, 0), false)
  def parse(j: Json): Either[JastError, Mask] = j match {
    case o: Json.Obj =>
      val ellp = if (o contains "ellipse") o("ellipse").to[Boolean].OUT     else default.ellipse
      val shap = if (o contains "shape")   o("shape").to[Array[Double]].OUT else default.shape.map(_.toDouble)
      val cutt = if (o contains "cut")     o("cut").to[Boolean].OUT         else default.cut
      if (shap.length != 4) return Left(JastError(f"Mask shape is given by four numbers, found ${shap.mkString(", ")}"))
      if (!shap.forall(x => x >= 0 && x.finite && sigfig(x).toInt == sigfig(x)))
        return Left(JastError(f"Mask coordinates must be integers, found ${shap.mkString(", ")}"))
      Right(new Mask(ellp, shap.map(x => sigfig(x).toInt), cutt))
    case _           => Left(JastError("Expected JSON object for Mask but did not get one"))
  }
}

case class Stimulus(
  device: String, id: String, port: String, channel: String, description: String,
  delay: Double, interval: Double, count: Int,
  high: Double, pulses: Int, pulseInterval: Double,
  shape: String,
  more: Option[Stimulus], custom: Json
) extends AsJson {
  private[this] def sigfig(x: Double) = math.rint(x*1e6) / 1e6 
  def json = Json ~
    ("device", device) ~ ("id", id) ~ ("port", port) ~ ("channel", channel) ~ ("description", description) ~
    ("delay", sigfig(delay)) ~ ("interval", sigfig(interval)) ~ ("count", count) ~
    ("high", sigfig(high)) ~ ("pulses", pulses) ~ ("pulse-interval", pulseInterval) ~
    ("shape", shape) ~? ("more", more) ~? ("custom", if (custom.isNull) None else Some(custom)) ~
    Json
}
object Stimulus extends FromJson[Stimulus] {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def default = new Stimulus("", "", "", "X", "default", 0, 0, 0, 0, 0, 0, "u", None, Json.Null)  
  def parse(j: Json): Either[JastError, Stimulus] = j match {
    case o: Json.Obj =>
      val devc = if (o contains "device")         o("device").to[String].OUT                    else default.device
      val iddd = if (o contains "id")             o("id").to[String].OUT                        else default.id
      val port = if (o contains "port")           o("port").to[String].OUT                      else default.port
      val chan = if (o contains "channel")        o("channel").to[String].OUT                   else default.channel
      val desc = if (o contains "description")    o("description").to[String].OUT               else default.description
      val dely = if (o contains "delay")          o("delay").to[Double].OUT.fn(sigfig)          else default.delay
      val intv = if (o contains "interval")       o("interval").to[Double].OUT.fn(sigfig)       else default.interval
      val cont = if (o contains "count")          o("count").to[Double].OUT.fn(sigfig)          else default.count
      val high = if (o contains "high")           o("high").to[Double].OUT.fn(sigfig)           else default.high
      val plsn = if (o contains "pulses")         o("pulses").to[Double].OUT .fn(sigfig)        else default.pulses
      val pint = if (o contains "pulse-interval") o("pulse-interval").to[Double].OUT.fn(sigfig) else default.pulseInterval
      val shap = if (o contains "shape")          o("shape").to[String].OUT                     else default.shape
      val more = if (o contains "more")           o("more").to(Stimulus).OUT.fn(x => Option(x)) else default.more
      val cust = o.get("custom") match { case Some(jj) => jj; case _ => default.custom }
      if (dely < 0 || !dely.finite || dely > 99999999) return Left(JastError(f"Bad value for delay: $dely"))
      if (intv < 0 || !intv.finite || intv > 99999999) return Left(JastError(f"Bad value for interval: $intv"))
      if (cont < 0 || cont.toInt != cont) return Left(JastError(f"Bad value for count: $cont"))
      if (high < 0 || !high.finite || high > 99999999) return Left(JastError(f"Bad value for stimulus high time: $high"))
      if (plsn < 0 || plsn.toInt != plsn) return Left(JastError(f"Bad value for pulse number: $plsn"))
      if (pint < 0 || !pint.finite || pint > 99999999) return Left(JastError(f"Bad value for pulse interval: $pint"))
      if (devc == "Ticklish") {
        if (!(Set("u", "i") contains shap)) return Left(JastError(f"Unknown stimulus shape: $shap"))
        if (chan.length != 1 || "ABCDEFGHIJKLMNOPQRSTUVWX".indexOf(chan) < 0)
          return Left(JastError(f"Bad channel: $chan"))
        if (cont > 0 && plsn > 0 && intv+1e-6 < sigfig(high + (plsn.toInt-1) * pint))
          return Left(JastError(f"Not enough time for stimuli in interval!"))
      }
      Right(new Stimulus(devc, iddd, port, chan, desc, dely, intv, cont.toInt, high, plsn.toInt, pint, shap, more, cust))
    case _           => Left(JastError("Expected JSON object for Stimulus Parameters but did not get one"))
  }
}

case class Coordinate(x: Double, y: Double) extends AsJson { def json = Json ~ ("x", x) ~ ("y" , y) ~ Json }
object Coordinate extends FromJson[Coordinate] {
  def parse(j: Json): Either[JastError, Coordinate] = {
      val x = j("x").double
      val y = j("y").double
      if (!x.finite || x < 0) return Left(JastError("No valid x coordinate"))
      if (!y.finite || y < 0) return Left(JastError("No valid y coordinate"))
      Right(new Coordinate(x, y))
  }
}

case class References(lowIntensity: Int, highIntensity: Int, coordinates: Array[Coordinate]) extends AsJson {
  def json = Json ~
    ("low-intensity", lowIntensity) ~ ("high-intensity", highIntensity) ~
    ("coordinates", Json(coordinates)) ~
    Json
}
object References extends FromJson[References] {
  val default = new References(0, 16000, new Array[Coordinate](0))
  def parse(j: Json): Either[JastError, References] = j match {
    case o: Json.Obj =>
      import ImplicitMwtSettingsJson.implicitMwtSettingsCoordinate
      val lowi = if (o contains "low-intensity")  o("low-intensity").to[Double].OUT          else default.lowIntensity
      val hiii = if (o contains "high-intensity") o("high-intensity").to[Double].OUT         else default.highIntensity
      val coor = if (o contains "coordinates")    o("coordinates").to[Array[Coordinate]].OUT else default.coordinates
      if (lowi < 0 || !lowi.finite) return Left(JastError(f"Low reference threshold must be finite and non-negative, was $lowi"))
      if (hiii < lowi || !hiii.finite) return(Left(JastError(f"High reference threshold must be finite and >= low $lowi, was $hiii")))
      Right(new References(math.rint(lowi).toInt, math.rint(hiii).toInt, coor))
    case _           => Left(JastError("Expected JSON object for References but did not get one"))
  }
}

case class CustomLabView(path: String, camera: String, warmup: Int, bitDepth: Int, autoStart: Double, aggregate: Boolean) extends AsJson { 
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  def json = Json ~ 
    ("path", path) ~  ("camera", camera) ~ ("warmup", warmup) ~
    ("bit-depth", bitDepth) ~ ("auto-start", sigfig(autoStart)) ~ ("aggregate", aggregate) ~
    Json
}
object CustomLabView extends FromJson[CustomLabView] {
  private[this] def sigfig(x: Double) = math.rint(x*1e3) / 1e3
  val default = new CustomLabView(".", "", 1, 0, 0)
  def parse(j: Json): Either[JastError, CustomLabView] = j match {
    case o: Json.Obj =>
      val path = if (o contains "path")       o("path").to[String].OUT       else default.path
      val camr = if (o contains "camera")     o("camera").to[String].OUT     else default.camera
      val warm = if (o contains "warmup")     o("warmup").to[Double].OUT     else default.warmup
      val bdep = if (o contains "bit-depth")  o("bit-depth").to[Double].OUT  else default.bitDepth
      val auts = if (o contains "auto-start") o("auto-start").to[Double].OUT else default.autoStart
      val aggr = if (o contains "aggregate")  o("aggregate").to[Boolean].OUT else default.aggregate
      if (!warm.finite || warm < 0 || warm.toInt != warm) return Left(JastError(f"Warmup must be a positive integer, was $warm"))
      if (!bdep.finite || bdep.toInt != bdep) return Left(JastError(f"Bit depth must be a non-negative integer, was $bdep"))
      if (!auts.finite || auts < 0) return Left(JastError(f"Auto-start time must be finite and non-negative, was $auts"))
      Right(new CustomLabView(path, camr, warm.toInt, bdep.toInt, sigfig(auts), aggr))
    case _           => Left(JastError("Expected JSON object for Custom LabView parameters but did not get one"))
  }
}

case class Settings(
  timestamp: Option[Either[java.time.LocalDateTime, java.time.Instant]],
  software: String, segmentation: Segmentation, output: Output,
  masks: Array[Mask], stimuli: Array[Stimulus], references: References,
  custom: Option[Either[Json, CustomLabView]]
) extends AsJson {
  private[this] def jstamp: Option[Json] = timestamp match {
    case Some(e) =>
      Some(e match {
        case Right(i) => Json(i)
        case Left(ldt) => Json(ldt)
      })
    case None => None
  }
  def json = (
    Json 
    ~? ("timestamp", jstamp) 
    ~ ("software", software) ~ ("segmentation", segmentation) ~ ("output", output)
    ~? ("masks", if (masks.length > 0) Some(Json(masks)) else None)
    ~? ("stimuli", if (stimuli.length > 0) Some(Json(stimuli)) else None)
    ~? ("references", if (references.coordinates.length > 0) Some(references) else None)
    ~? ("custom", custom.map(_ match { case Right(clv) => Json(clv); case Left(j) => j }))
    ~ Json
  )
  override def toString = PrettyJson(json)
}
object Settings extends FromJson[Settings] {
  val default = new Settings(None, "", Segmentation.default, Output.default, new Array[Mask](0), new Array[Stimulus](0), References.default, None)
  def parse(j: Json): Either[JastError, Settings] = j match {
    case o: Json.Obj =>
      import ImplicitMwtSettingsJson._
      val time = if (o contains "timestamp")    o("timestamp").to[Either[java.time.LocalDateTime, java.time.Instant]].OUT else null
      val soft = if (o contains "software")     o("software").to[String].OUT                                              else default.software
      val segm = if (o contains "segmentation") o("segmentation").to[Segmentation].OUT                                    else default.segmentation
      val outp = if (o contains "output")       o("output").to[Output].OUT                                                else default.output
      val mask = if (o contains "masks")        o("masks").to[Array[Mask]].OUT                                            else default.masks
      val stim = if (o contains "stimuli")      o("stimuli").to[Array[Stimulus]].OUT                                      else default.stimuli
      val refs = if (o contains "references")   o("references").to[References].OUT                                        else default.references
      val cust = o("custom") match {
        case j: Json => Some(j.to[CustomLabView] match { case Right(clv) => Right(clv); case _=> Left(j) })
        case _       => None
      }
      Right(new Settings(Option(time), soft, segm, outp, mask, stim, refs, cust))
    case Json.Null   => Right(default)
    case _           => Left(JastError("Expected JSON object for MWT settings but did not find one"))
  }
}
