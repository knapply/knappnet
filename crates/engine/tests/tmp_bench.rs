use knappnet::metrics::partition::*;
use knappnet::Graph;
fn load(name: &str) -> Graph {
    let txt = std::fs::read_to_string(format!("../../validation/fixtures/{name}/edges.csv")).unwrap();
    let meta = std::fs::read_to_string(format!("../../validation/fixtures/{name}/meta.json")).unwrap();
    let n: usize = meta.split("\"n\": ").nth(1).unwrap().split(',').next().unwrap().trim().parse().unwrap();
    let directed = meta.contains("\"directed\": true");
    let (mut s, mut t, mut w) = (vec![], vec![], vec![]);
    for line in txt.lines().skip(1) {
        let f: Vec<&str> = line.split(',').collect();
        s.push(f[0].parse::<i64>().unwrap()); t.push(f[1].parse::<i64>().unwrap());
        w.push(f.get(2).and_then(|x| x.parse::<f64>().ok()).unwrap_or(1.0));
    }
    Graph::from_edges(n, &s, &t, Some(&w), directed).unwrap()
}
#[test]
fn bench() {
    for name in std::env::var("FX").unwrap().split(",") {
        let g = load(name);
        let t0 = std::time::Instant::now();
        let r = optimal_modularity(&g, 1.0, false);
        let wt = walktrap(&g.clone(), 4, false, 1).ok().map(|w| w.modularity);
        println!("{name}: {:?} {:?} walktrap {:?}", r.map(|r| r.modularity), t0.elapsed(), wt);
    }
}
