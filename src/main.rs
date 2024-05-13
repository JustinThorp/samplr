use std::f64::consts::PI;

use rand::Rng;


trait Distribution {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> f64;
    fn pdf(&self,x:f64) -> f64;
    fn lnpdf(&self,x: f64) -> f64;
}

fn main() {
    // let mut rng = rand::thread_rng();
     let x: Normal = Normal::new(100.0,1.0);
    // let y: Vec<f64> = (1..10000).into_iter().map(|_| x.sample(&mut rng)).collect();
    // let mu: f64 = y.iter().sum::<f64>() / y.len() as f64;
    //println!("{:?}",x.lnpdf(10.0));
    merto_hast();
}


fn merto_hast () {
    let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
    let mut x: Vec<f64> = Vec::new();
    x.push(0.0);
    let p = Normal::new(100.0,1.0);
    let mut j = 0;
    for i in 1..2000 {
        let g = Normal::new(x[i-1],1.0);
        let t = g.sample(&mut rng);
        //let a = p.pdf(t) / p.pdf(x[i-1]);
        let a = (p.lnpdf(t) - p.lnpdf(x[i-1])).exp();
        let u: f64 = rng.gen();
        //println!("{a};{u}");
        if u <= a {
            j += 1;
            x.push(t);
        } else {
            x.push(x[i-1])
        }
    }

    println!("{:?}",j as f64 / 2000.0);
}


struct Normal {
    mu: f64,    
    sigma: f64
}


impl Normal {
    fn new(mu: f64,sigma: f64) -> Normal{
        Normal { mu, sigma}
    }

}

impl Distribution for  Normal {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> f64 
    {
        let u: f64 = rng.gen();
        let v: f64 = rng.gen();
        let x: f64 =(-2.0*u.ln()).sqrt() * (2.0*PI*v).cos();
        return x*self.sigma + self.mu;
    }


    fn pdf(&self,x: f64) -> f64 {
        //(1.0 / (self.sigma*(PI*2.0).sqrt())) * (-0.5*((x - self.mu) / self.sigma).powi(2)).exp()
        (-0.5*((x - self.mu) / self.sigma).powi(2)).exp()
    }

    fn lnpdf(&self,x: f64) -> f64 {
        -1.0*(self.sigma).ln() - 0.5*(2.0*PI).ln() - 0.5*((x - self.mu) / self.sigma).powi(2)
    }
    
    
}

