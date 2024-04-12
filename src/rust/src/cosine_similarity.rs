pub fn cosine_similarity(x: Vec<f64>,
          y: Vec<f64>
         ) -> f64 {

    // dot product
    let p: Vec<f64> = x.clone().into_iter().zip(y.clone()).map(|(x, y)| x * y).collect();
    let p: f64 = p.iter().sum::<f64>() as f64;
    // x^2 > sum > sqrt
    let xd: Vec<f64> = x.clone().into_iter().map(|x| x.powf(2.0)).collect();
    let xd: f64 =  xd.iter().sum::<f64>().sqrt() as f64;
    // y^2 > sum > sqrt
    let yd: Vec<f64> = y.clone().into_iter().map(|x| x.powf(2.0)).collect();
    let yd: f64 =  yd.iter().sum::<f64>().sqrt() as f64;

    // cosine similarity
    return p / (xd * yd);
}
