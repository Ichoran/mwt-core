extern crate clap;
extern crate cgmath;

use std::fmt::Display;
use std::ops::*;
use std::str::FromStr;
use clap::{Arg, ArgMatches, App};
use cgmath::{Vector2};

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


#[derive(Debug)]
struct Worm {
    ampl: f64,
    wave: f64,
    cycles: f64,
    aspect: f64,
    origin: Vector2<f32>,
    step: Vector2<f32>,
    shape: Box<fn(f64) -> f64>,
}

fn widely(position: f64) -> f64 {
    if      position < 0.1 { position }
    else if position < 0.3 { 0.1 + (position-0.1)*4.5 }
    else if position < 0.6 { 1.0 }
    else                   { (1.0 - position)*2.5 }
}

impl Worm {
    fn new(amp: f64, wav: f64, cyc: f64, asp: f64, ori: Vector2<f32>, ste: Vector2<f32>, sha: Box<fn(f64) -> f64>) -> Worm {
        Worm{ ampl: amp, wave: wav, cycles: cyc, aspect: asp, origin: ori, step: ste, shape: sha }
    }
    fn default() -> Worm {
        Worm { 
            ampl: 4.0, wave: 20.0, cycles: 1.5, aspect: 0.1,
            origin: Vector2::new(16f32, 16f32),
            step:   Vector2::new( 1f32,  1f32),
            shape: Box::new(widely)
        }
    }
    fn parse(string: &str) -> Ret<Worm> {
        let mut worm = Worm::default();
        let mut tok = string.split(':');
        let ampl = pf64(tok.next().ok_or("Worm has no amplitude")?, "worm amplitude")?;
        let wave = pf64(tok.next().ok_or("Worm has no wavelength")?, "worm wavelength")?;
        let cycl = pf64(tok.next().ok_or("Worm has no cycle length")?, "worm cycles")?;
        let aspc = pf64(tok.next().ok_or("Worm has no aspect ratio")?, "worm aspect")?;
        let orig = pvf32(tok.next().ok_or("Worm has no origin vector")?, "worm origin")?;
        let step = pvf32(tok.next().ok_or("Worm has no movement axis")?, "worm step")?;
        Ok(Worm::new(ampl, wave, cycl, aspc, orig, step, Box::new(widely)))
    }
    fn at_time<'a>(&'a self, t: f64) -> WormIter<'a> {
        WormIter{ time: t, axial: -1f64, radial: 0f64, worm: self }
    }
}

#[derive(Debug)]
struct WormIter<'a> {
    time: f64,
    axial: f64,
    radial: f64,
    worm: &'a Worm
}

impl<'a> Iterator for WormIter<'a> {
    type Item = Vector2<f32>;
    fn next(&mut self) -> Option<Self::Item> { None }
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

impl Image {
    fn new(w: u32, h: u32, bg: u8) -> Image {
        let size = (w*h) as usize;
        Image { width: w, height: h, data: vec![bg; size] }
    }

    fn imprint(&mut self, fg: u8, worm: Worm, t: f64) -> &mut Image {
        self[Vector2{ x: 1, y: 1}] = fg;
        self[Vector2{ x: 2, y: 1}] = fg;
        self
    }
}

fn main() {
    let args = App::new("Fake Worm Generator")
        .version("0.1")
        .author("Rex A. Kerr (ichoran@gmail.com)")
        .about("Generates a series of raw image files with worm-like shapes embedded")
        .arg(Arg::with_name("width") .short("w").long("width") .takes_value(true))
        .arg(Arg::with_name("height").short("h").long("height").takes_value(true))
        .arg(Arg::with_name("frames").short("f").long("frames").takes_value(true))
        .arg(Arg::with_name("bg")               .long("bg")    .takes_value(true))
        .arg(Arg::with_name("fg")               .long("fg")    .takes_value(true))
        .arg(Arg::with_name("noise")            .long("noise") .takes_value(true))
        .arg(Arg::with_name("out")   .short("o").long("output").takes_value(true))
        .arg(Arg::with_name("WORMS").index(1).max_values(1000).help(
            "Specify one or more animal positions in format ampl:wave:len:aspect:px^py:vx^vy"
        ))
        .get_matches();
    let params = match Parameters::parse(&args) {
        Err(e) => { println!("Error in arguments:\n{}", e); std::process::exit(1) }
        Ok(p)  => p
    };
    println!("{:?}", params);
}
