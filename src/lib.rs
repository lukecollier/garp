use anyhow::Result;
use geo::geometry::Point;
use geo::prelude::*;
use ordered_float::OrderedFloat;
use rand::Rng;
use uom::si::f32::*;
use uom::si::length::meter;

struct Vehicle {}
#[derive(Clone, PartialEq, Copy, Debug, PartialOrd)]
struct Location {
    lat: f32,
    lng: f32,
}

impl Location {
    fn new() -> Location {
        Location {
            lat: 0_f32,
            lng: 0_f32,
        }
    }
    fn random() -> Location {
        let mut rand = rand::thread_rng();
        let lat = rand.gen_range(-180_f32..180_f32);
        let lng = rand.gen_range(-90_f32..90_f32);
        Location { lat, lng }
    }
}
struct Resource {
    can_use_vehicle: Vec<Vehicle>,
    start_location: Vec<Vehicle>,
    end_location: Vec<Vehicle>,
}
// introduce capex https://onlinelibrary.wiley.com/doi/full/10.1002/net.22028
// introduces VRP and a lot of data about cluster first route second, and route first cluster
// second https://www.mdpi.com/2079-9292/10/24/3147
struct Segment {
    from: Location,
    to: Location,
}

struct TravelingSalesmanProblem {
    visits: Vec<Location>,
    distances: Vec<Segment>,
    vehicles: Vec<Vehicle>,
}

struct VehicleRoutingProblem {
    visits: Vec<Location>,
    distances: Vec<Segment>,
    vehicles: Vec<Vehicle>,
}
// todo: Could do this as a binary string as described here https://www.mdpi.com/2079-9292/10/24/3147#B5-electronics-10-03147
#[derive(Clone)]
struct Solution {
    route: Vec<Location>,
}

// polymorphic way to pass in a cost
trait Costs {
    fn between(&self, one: &Location, two: &Location) -> Option<Length>;
}

struct Distances {
    distance_matrix: Vec<Vec<Length>>,
    from: Vec<Location>,
    to: Vec<Location>,
}

impl Costs for Distances {
    fn between(&self, loc_one: &Location, loc_two: &Location) -> Option<Length> {
        let from_pos = self.from.iter().position(|loc| loc_one == loc);
        let to_pos = self.to.iter().position(|loc| loc_two == loc);
        from_pos.and_then(|from| to_pos.map(|to| self.distance_matrix[from][to]))
    }
}

struct HaversineDistances {}
impl HaversineDistances {
    fn new() -> HaversineDistances {
        HaversineDistances {}
    }
}

impl Costs for HaversineDistances {
    fn between(&self, loc_one: &Location, loc_two: &Location) -> Option<Length> {
        let point_one = Point::new(loc_one.lat, loc_one.lng);
        let point_two = Point::new(loc_two.lat, loc_two.lng);
        Some(Length::new::<meter>(
            point_one.haversine_distance(&point_two),
        ))
    }
}

//Roulette wheel selection (RWS)—chances of an individual being chosen are proportional to its fitness value; thus, selection may be imagined as a spinning roulette, where each individual takes an amount of space on the roulette wheel according to its fitness.
//Elitism selection (ES)—a certain percentage of the population, ordered by fitness, is always transferred to the next population. In that scenario, the algorithm makes sure that best so far known solutions would not be lost in the process of selection.
//Rank selection (RS)—similar to RWS, but each individual solution’s space on the roulette wheel is not proportional to its fitness, but to its rank in the list of all individuals, ordered by fitness.
//Stochastic universal sampling selection (SUSS)—instead of spinning the wheel of the roulette for a certain amount of times, spin it once. If selecting n individuals, there must exist n spaces on the wheel, and the chosen individual is copied n times to the next generation.
//Tournament selection (TS)—as many times as required, choose two individuals randomly, and let the more fit one be chosen.
fn fitness_function(solution: &Solution, distances: &Distances) -> usize {
    todo!()
}

// Technically undermines calling this GA, but hey ho here we go
fn cost_function<C: Costs>(solution: &Solution, costs: &C) -> f32 {
    let total_distance: Length = solution
        .route
        .windows(2)
        .into_iter()
        .map(|loc| (costs.between(&loc[0], &loc[1])).unwrap())
        .sum();
    total_distance.get::<meter>()
}
// plan is apply this first, the returned population is left (the alive) and the right (killed)
// take the right populatio and apply it to the selection tournament, then we're sorted.
// I think we'll do elitism + tournament selection
// and then we also
// todo: make these selection chains
fn selection_elitism(solutions: &Vec<Solution>) -> (Vec<Solution>, Vec<Solution>) {
    let haversine_distances_costs = HaversineDistances::new();
    let mut ordered_solutions: Vec<(&Solution, OrderedFloat<f32>)> = solutions
        .into_iter()
        .map(|solution| {
            (
                solution,
                OrderedFloat::from(cost_function(&solution, &haversine_distances_costs)),
            )
        })
        .collect::<Vec<_>>();
    ordered_solutions.sort_by(|(_, first), (_, second)| first.cmp(second));

    let mid = (solutions.len() / 100) * 10;
    let (left, right) = ordered_solutions
        .into_iter()
        .map(|(solution, _)| solution.to_owned())
        .collect::<Vec<_>>()
        .split_at(mid);
    todo!()
}
fn selection_tournament(population: Vec<Solution>) -> (Vec<Solution>, Vec<Solution>) {
    let haversine_distances_costs = HaversineDistances::new();
    todo!()
}
fn selection(population: Vec<Solution>) -> Vec<Solution> {
    todo!()
}

// Take both parents and randomly choose two crossover points, the same for both of them.
// Copy the integers in between the crossover points, from male parent to child, keeping them at the same positions.
// Take the female parent, and starting from the gene after the second crossover point, iterate through all genes. If the end is met, start from the beginning.
// Take the child and starting from the gene after the second crossover point, copy the female parent gene that is considered in the current iteration, only if it is not present yet in the child’s chromosome. If the end is met, start from the beginning.
// The operation is finished if all empty spaces in the child chromosome are filled.
// Optionally, swap the roles of female and male parents, and repeat the whole process to produce a second offspring.
fn order_crossover(parent_one: &Vec<Location>, parent_two: &Vec<Location>) -> Vec<Location> {
    let mut rand = rand::thread_rng();
    let decider = rand.gen::<bool>();
    let (male, female) = if decider {
        (parent_one, parent_two)
    } else {
        (parent_two, parent_one)
    };
    assert!(&female.len() == &male.len());
    let start_at = rand.gen_range(0..male.len());
    let end_at = rand.gen_range(start_at..male.len());

    let first_second = &female[start_at..female.len()];
    let female_second = &female[0..end_at];
    let mut female_chromosomes = [first_second, female_second].concat();

    let mut child = vec![Location::new(); male.len()];
    child.splice(start_at..end_at, male[start_at..end_at].to_owned());

    for idx in end_at..child.len() {
        // prune the list
        // future: replace with https://doc.rust-lang.org/std/vec/struct.Vec.html#method.drain_filter
        let mut i = 0;
        while i < female_chromosomes.len() {
            if child.contains(&female_chromosomes[i]) {
                female_chromosomes.remove(i);
            } else {
                i += 1;
            }
        }
        let _ = std::mem::replace(&mut child[idx], female_chromosomes.pop().unwrap());
    }

    for idx in 0..start_at {
        // prune the list
        // future: replace with https://doc.rust-lang.org/std/vec/struct.Vec.html#method.drain_filter
        let mut i = 0;
        while i < female_chromosomes.len() {
            if child.contains(&female_chromosomes[i]) {
                female_chromosomes.remove(i);
            } else {
                i += 1;
            }
        }
        let _ = std::mem::replace(&mut child[idx], female_chromosomes.pop().unwrap());
    }
    assert!(child.len() == male.len());
    assert!(female_chromosomes.len() == 0);
    child.into()
}

// takes in a population and produces a new population with children
fn reproduction(population: Vec<Solution>) -> Vec<Solution> {
    todo!()
}

// Take the population given and mutate it via swapping routes from each solution
// Could we mutate using the distance matrix as well? I.e take two points and choose a different
// route (if available)
fn mutation(population: Vec<Solution>) -> Vec<Solution> {
    todo!("here we mutate the population to avoid local minima")
}

#[cfg(test)]
mod tests {
    use assert_unordered::assert_eq_unordered;

    use super::*;
    use rand::prelude::*;

    #[test]
    fn order_crossover_children_contains_parent_paths() {
        let mut path = vec![Location::new(); 10];
        path.fill_with(|| Location::random());
        let mut rand = rand::thread_rng();
        let mut parent_one = path.clone();
        parent_one.shuffle(&mut rand);
        let mut parent_two = path.clone();
        parent_two.shuffle(&mut rand);
        let child = order_crossover(&parent_two, &parent_one);
        assert_eq_unordered!(child, path);
    }

    #[test]
    fn order_crossover_children_is_consistent() {
        let mut path = vec![Location::new(); 10];
        path.fill_with(|| Location::random());
        let mut rand = rand::thread_rng();
        let mut parent_one = path.clone();
        parent_one.shuffle(&mut rand);
        let mut parent_two = path.clone();
        parent_two.shuffle(&mut rand);
        let child = order_crossover(&parent_two, &parent_one);
        assert_eq_unordered!(child, path);
    }
}
