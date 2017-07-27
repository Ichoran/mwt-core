extern crate clap;
extern crate cgmath;

use std::ops::*;
use clap::{Arg, ArgMatches, App};
use cgmath::{Vector2};

struct Worm {
    ampl: f64,
    wave: f64,
    cycles: f64,
    origin: Vector2<f32>,
    step: Vector2<f32>,
    shape: Box<fn(f64) -> f64>,
}

struct Parameters {
    width: u32,
    height: u32,
    frames: u32,
    background: u8,
    foreground: u8,
    noise: u8,
    prefix: String,
    defaults: Worm,
    worms: Vec<Worm>
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

impl Image {
    fn new(w: u32, h: u32, bg: u8) -> Image {
        let size = (w*h) as usize;
        Image { width: w, height: h, data: vec![bg; size] }
    }
}

fn widely(position: f64) -> f64 {
    if      position < 0.1 { position }
    else if position < 0.3 { 0.1 + (position-0.1)*4.5 }
    else if position < 0.6 { 1.0 }
    else                   { (1.0 - position)*2.5 }
}

fn parse_parameters(args: &ArgMatches)-> Parameters {
    let default_worm = 
        Worm { 
            ampl: 4.0, wave: 20.0, cycles: 1.5,
            origin: Vector2::new(16f32, 16f32),
            step:   Vector2::new( 1f32,  1f32),
            shape: Box::new(widely)
        };
    let mut default = Parameters {
        width: 64, height: 64, frames: 16,
        background: 240, foreground: 190, noise: 10,
        prefix: "".to_string(),
        defaults: default_worm,
        worms: Vec::new()
    };
    default
}

fn main() {
    println!("Hello, world!");
    let args = App::new("Fake Worm Generator")
        .version("0.1")
        .author("Rex A. Kerr (ichoran@gmail.com)")
        .about("Generates a series of raw image files with worm-like shapes embedded")
        .arg(Arg::with_name("width").short("w").long("width").takes_value(true))
        .arg(Arg::with_name("height").short("h").long("height").takes_value(true))
        .arg(Arg::with_name("frames").short("f").long("frames").takes_value(true))
        .arg(Arg::with_name("bg").long("bg").takes_value(true))
        .arg(Arg::with_name("fg").long("fg").takes_value(true))
        .arg(Arg::with_name("noise").long("noise").takes_value(true))
        .arg(Arg::with_name("out").short("o").long("output").takes_value(true))
        .arg(Arg::from_usage("<WORMS>").max_values(1000000))
        .get_matches();
    let params = parse_parameters(&args);
    println!("Goodbye, world!");
}
