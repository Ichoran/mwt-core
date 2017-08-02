extern crate clap;
extern crate cgmath;
extern crate rand;

use std::fmt::Display;
use std::io::Write;
use std::ops::*;
use std::str::FromStr;
use clap::{Arg, ArgMatches, App};
use cgmath::{Vector2, InnerSpace};

type Ret<T> = Result<T, String>;

fn px<T, E>(s: &str, msg: &str) -> Ret<T>
where T: FromStr<Err=E>, E: Display {
    s.parse::<T>().map_err(|e| {
        if msg.len() > 0 { format!("{} in {}", e, msg) }
        else             { format!("{}", e) }
    })
}

fn pu8(s: &str, msg: &str) -> Ret<u8>   { px(s, msg) }

fn pu32(s: &str, msg: &str) -> Ret<u32> { px(s, msg) }

fn pf64(s: &str, msg: &str) -> Ret<f64> { px(s, msg) }

fn pvf32(s: &str, msg: &str) -> Ret<Vector2<f32>> {
    let mut tok = s.split('^');
    let x = pf64(tok.next().ok_or("No X component")?, msg)?;
    let y = pf64(tok.next().ok_or("No Y component")?, msg)?;
    Ok(Vector2::new(x as f32, y as f32))
}


#[derive(Clone, Copy, Debug)]
struct Worm {
    ampl: f64,
    wave: f64,
    cycles: f64,
    aspect: f64,
    phase: f64,
    origin: Vector2<f32>,
    step: Vector2<f32>,
}

fn widely(position: f64) -> f64 {
    if      position < 0.2 { position*3.0 }
    else if position < 0.5 { 0.6 + (position - 0.2)*1.333333333333 }
    else if position < 0.7 { 1.0 }
    else                   { (1.0 - position)*3.333333333333333 }
}

impl Worm {
    fn new(amp: f64, wav: f64, cyc: f64, asp: f64, pha: f64, ori: Vector2<f32>, ste: Vector2<f32>) -> Worm {
        Worm{ ampl: amp, wave: wav, cycles: cyc, aspect: asp, origin: ori, phase: pha, step: ste }
    }
    fn default() -> Worm {
        Worm { 
            ampl: 4.0, wave: 20.0, cycles: 1.5, aspect: 0.1, phase: 0.0,
            origin: Vector2::new(16f32, 16f32),
            step:   Vector2::new( 1f32,  1f32)
        }
    }
    fn parse(string: &str) -> Ret<Worm> {
        let mut tok = string.split(':');
        let ampl = pf64(tok.next().ok_or("Worm has no amplitude")?, "worm amplitude")?;
        let wave = pf64(tok.next().ok_or("Worm has no wavelength")?, "worm wavelength")?;
        let cycl = pf64(tok.next().ok_or("Worm has no cycle length")?, "worm cycles")?;
        let aspc = pf64(tok.next().ok_or("Worm has no aspect ratio")?, "worm aspect")?;
        let phas = pf64(tok.next().ok_or("Worm has no phase")?, "worm phase")?;
        let orig = pvf32(tok.next().ok_or("Worm has no origin vector")?, "worm origin")?;
        let step = pvf32(tok.next().ok_or("Worm has no movement axis")?, "worm step")?;
        let ans = Worm::new(ampl, wave, cycl, aspc, phas, orig, step);
        Ok(ans)
    }
    fn at_time<'a>(&'a self, t: f64) -> WormIter<'a> {
        fn zero() -> Vector2<f32> { Vector2::new(0f32, 0f32) }
        let mut ans = WormIter{ time: t, axial: 0f64, radial: 0f64, spine: zero(), orth: zero(), worm: self };
        ans.locate();
        ans
    }
}

#[derive(Debug)]
struct WormIter<'a> {
    time: f64,
    axial: f64,
    radial: f64,
    spine: Vector2<f32>,
    orth: Vector2<f32>,
    worm: &'a Worm
}

impl<'a> WormIter<'a> {
    fn locate(&mut self) {
        let L = self.worm.step.magnitude() as f64;
        let u = self.worm.step.normalize();
        let v = Vector2::new(-u.y, u.x);
        let x = (self.axial + self.time*L) as f32;
        let theta = 2f32 * std::f32::consts::PI * x / (self.worm.wave as f32) + (self.worm.phase as f32);
        let y = (self.worm.ampl as f32) * theta.sin();
        let wiggle = self.worm.wave/(2.0 * std::f64::consts::PI * self.worm.ampl);
        self.spine = self.worm.origin + u*x + v*y;
        self.orth = Vector2::new(-theta.cos(), wiggle as f32).normalize();
    }
}

impl<'a> Iterator for WormIter<'a> {
    type Item = Vector2<f32>;
    fn next(&mut self) -> Option<Self::Item> {
        let L = self.worm.cycles * self.worm.wave;
        if self.axial > 0.5*L { return None; }
        let w = 0.5*L * self.worm.aspect * widely(self.axial/L + 0.5);
        let pt = self.spine + (self.orth * (self.radial as f32));
        if self.radial <= 0.0 { 
            self.radial -= 0.1;
            if self.radial < -w { self.radial = 0.1; }
        }
        else { self.radial += 0.1; }
        if self.radial > w {
            self.radial = 0.0;
            if self.axial <= 0.0 {
                self.axial -= 0.1;
                if self.axial < -0.5*L {
                    self.axial = 0.1;
                }
            }
            else { self.axial += 0.1; }
            self.locate();
        }
        Some(pt)
    }
}

#[derive(Debug)]
struct Parameters {
    width: u32,
    height: u32,
    frames: u32,
    background: u8,
    foreground: u8,
    noise: u8,
    prefix: String,
    worms: Vec<Worm>
}

impl Parameters {
    fn parse(args: &ArgMatches)-> Ret<Parameters> {
        let mut ans = Parameters {
            width: 64, height: 64, frames: 16,
            background: 240, foreground: 190, noise: 10,
            prefix: "".to_string(),
            worms: Vec::new()
        };
        if let Some(w) = args.value_of("width")  { ans.width  = pu32(w, "width")? }
        if let Some(h) = args.value_of("height") { ans.height = pu32(h, "height")? }
        if let Some(f) = args.value_of("frames") { ans.frames = pu32(f, "frames")? }
        if let Some(bg)= args.value_of("bg")     { ans.background = pu8(bg, "bg")? }
        if let Some(fg)= args.value_of("fg")     { ans.foreground = pu8(fg, "fg")? }
        if let Some(n) = args.value_of("noise")  { ans.noise  = pu8(n, "noise")? }
        if let Some(s) = args.value_of("prefix") { ans.prefix = s.to_string() }
        if let Some(ws)= args.values_of("WORMS") {
            for (i,w) in ws.enumerate() {
                ans.worms.push(Worm::parse(w).map_err(|e| format!("Error in worm {}\n{}", i+1, e))?);
            }
        }
        Ok(ans)
    }    
}

struct Image {
    width: u32,
    height: u32,
    data: Vec<u8>
}

impl Index<Vector2<u32>> for Image {
    type Output = u8;
    fn index(&self, index: Vector2<u32>) -> &u8 {
        let x = index.x % self.width;
        let y = index.y % self.height;
        let data: &[u8] = self.data.as_ref();
        &data[(x + y*self.width) as usize]
    }
}

impl IndexMut<Vector2<u32>> for Image {
    fn index_mut(&mut self, index: Vector2<u32>) -> &mut u8 {
        let x = index.x % self.width;
        let y = index.y % self.height;
        let mut data: &mut [u8] = self.data.as_mut();
        &mut data[(x + y*self.width)as usize]
    }
}

impl<'a> From<&'a ImageF32> for Image {
    fn from(that: &ImageF32) -> Self {
        let data: Vec<u8> = that.data.iter().map(|x| if *x > 255f32 { 255u8 } else if *x > 0f32 { *x as u8 } else { 0u8 }).collect();
        Image{ width: that.width, height: that.height, data: data }
    }
}

impl Image {
    fn new(w: u32, h: u32, bg: u8) -> Image {
        let size = (w*h) as usize;
        Image { width: w, height: h, data: vec![bg; size] }
    }
    fn ascii(&self) -> String {
        let mut s = String::new();
        let mut hi = self[Vector2::new(0, 0)];
        let mut lo = self[Vector2::new(0, 0)];
        for i in self.data.iter() {
            if *i > hi { hi = *i; }
            if *i < lo { lo = *i; }
        }
        let table = vec!["  ", " .", "..", ".=", "==", "=X", "XX", "X#", "##", "@@"];
        let l = lo as f64;
        let h = hi as f64;
        for y in 0..self.height {
            for x in 0..self.width {
                let i = self[Vector2::new(x, y)] as f64;
                let v = (10.0*(i - l) / (1.0 + h - l)) as usize;
                s.push_str(table[9 - v]);
            }
            s.push('\n');
        }
        s
    }
}

struct ImageF32 {
    width: u32,
    height: u32,
    background: f32,
    foreground: f32,
    data: Vec<f32>
}

impl Index<Vector2<u32>> for ImageF32 {
    type Output = f32;
    fn index(&self, index: Vector2<u32>) -> &f32 {
        let x = index.x % self.width;
        let y = index.y % self.height;
        let data: &[f32] = self.data.as_ref();
        &data[(x + y*self.width) as usize]
    }
}

impl IndexMut<Vector2<u32>> for ImageF32 {
    fn index_mut(&mut self, index: Vector2<u32>) -> &mut f32 {
        let x = index.x % self.width;
        let y = index.y % self.height;
        let mut data: &mut [f32] = self.data.as_mut();
        &mut data[(x + y*self.width)as usize]
    }
}

impl ImageF32 {
    fn new(w: u32, h: u32, bg: f32, fg: f32) -> ImageF32 {
        let size = (w*h) as usize;
        ImageF32 { width: w, height: h, background: bg, foreground: fg, data: vec![bg; size] }
    }

    fn blank(&mut self) {
        for i in self.data.iter_mut() { *i = self.background }
    }

    fn noisify(&mut self, sd: f64) {
        for i in self.data.iter_mut() {
            let rand::distributions::normal::StandardNormal(x) = rand::random();
            *i += (x * sd) as f32;
        }
    }

    fn imprint(&mut self, worm: Worm, t: f64) -> &mut ImageF32 {
        let delta = (self.foreground - self.background)*0.02;
        for pt in worm.at_time(t) {
            let ij = Vector2::new(pt.x.round() as u32, pt.y.round() as u32);
            let i = self[ij] + delta;
            self[ij] = { if (self.foreground - i)*delta < 0.0 { self.foreground } else { i } };
        }
        self
    }
}

fn main() {
    let args = App::new("Fake Worm Generator")
        .version("0.1")
        .author("Rex A. Kerr (ichoran@gmail.com)")
        .about("Generates a series of raw image files with worm-like shapes embedded")
        .arg(Arg::with_name("width") .short("w").long("width") .takes_value(true).help("Width of image in pixels"))
        .arg(Arg::with_name("height").short("h").long("height").takes_value(true).help("Height of image in pixels"))
        .arg(Arg::with_name("frames").short("f").long("frames").takes_value(true).help("Number of images"))
        .arg(Arg::with_name("bg")               .long("bg")    .takes_value(true).help("Background intensity (0-255)"))
        .arg(Arg::with_name("fg")               .long("fg")    .takes_value(true).help("Foreground intensity (0-255)"))
        .arg(Arg::with_name("noise")            .long("noise") .takes_value(true).help("Noise (standard deviation)"))
        .arg(Arg::with_name("out")   .short("o").long("output").takes_value(true).help("Base filename (no extension)"))
        .arg(Arg::with_name("WORMS").index(1).max_values(1000).help(
            "Specify one or more animal positions in format ampl:wave:len:aspect:phase:px^py:vx^vy"
        ))
        .get_matches();
    let params = match Parameters::parse(&args) {
        Err(e) => { println!("Error in arguments:\n{}", e); std::process::exit(1) }
        Ok(p)  => p
    };
    println!("{:?}", params);
    let mut fim = ImageF32::new(params.width, params.height, params.background as f32, params.foreground as f32);
    let mut imgs = Vec::with_capacity(params.frames as usize);
    for n in 0..params.frames {
        fim.blank();
        for w in &params.worms {
            fim.imprint(*w, n as f64);
        }
        fim.noisify(params.noise as f64);
        imgs.push(Image::from(&fim));
    }
    if params.frames > 0 && params.width < 40 {
        println!("{}", imgs[0].ascii());
    }
    let digits = (format!("{}", if params.frames > 0 { params.frames-1 } else { 0 })).len();
    let name = |n: u32| {
        let mut s = String::new();
        s.push_str(params.prefix.as_str());
        let num = format!("{}", if n >= params.frames { params.frames-1 } else { n });
        for _ in 0..(digits - num.len()) { s.push('0'); }
        s.push_str(num.as_str());
        s.push_str(".raw");
        s
    };
    for n in 0..params.frames {
        let mut f = std::fs::File::create(name(n)).expect("Could not create file");
        f.write_all(imgs[n as usize].data.as_slice()).expect("Error writing");
    }
}
