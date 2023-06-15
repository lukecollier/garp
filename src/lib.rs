use anyhow::Result;
use geo::geometry::Point;
use geo::prelude::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use percentage::Percentage;
use percentage::PercentageInteger;
use rand::seq::SliceRandom;
use rand::Rng;
use uom::si::f32::*;
use uom::si::length::meter;
use uom::Conversion;

struct Vehicle {}
#[derive(Clone, PartialEq, Copy, Debug, PartialOrd)]
pub struct Location {
    lat: f32,
    lng: f32,
}

impl Location {
    pub fn new() -> Location {
        Location {
            lat: 0_f32,
            lng: 0_f32,
        }
    }
    pub fn random() -> Location {
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
pub struct Solution {
    route: Vec<Location>,
}

impl Solution {
    pub fn genetic_algorithm(
        population: Vec<Solution>,
        generations_limit: usize,
        success_predicate: &dyn Fn(f32) -> bool,
    ) -> Vec<Solution> {
        let mut current_generation = Vec::new();
        for _ in 0..generations_limit {
            current_generation = next_generation(&population);
        }
        current_generation
    }
}

// costs r like this as costs will want to encapsulate state and structs let us do that bby
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

// probably pass this in to the functions as a dynamic function
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
    let solutions = ordered_solutions
        .into_iter()
        .map(|(solution, _)| solution.clone())
        .collect::<Vec<_>>();
    let (selected, not_selected) = solutions.split_at(mid);
    (selected.to_vec(), not_selected.to_vec())
}
fn selection_tournament(population: &Vec<Solution>) -> (Vec<Solution>, Vec<Solution>) {
    let haversine_distances_costs = HaversineDistances::new();
    todo!()
}
// do some DI for the cost fn, probably wants to be &dyn Fn(&solution) -> f32
fn selection(population: &Vec<Solution>) -> (Vec<Solution>, Vec<Solution>) {
    let (mut elitists, mut not_elitists) = selection_elitism(&population);
    let (mut winners, mut losers) = selection_tournament(&not_elitists);
    losers.append(&mut not_elitists);
    winners.append(&mut elitists);
    (winners, losers)
}

// Take both parents and randomly choose two crossover points, the same for both of them.
// Copy the integers in between the crossover points, from male parent to child, keeping them at the same positions.
// Take the female parent, and starting from the gene after the second crossover point, iterate through all genes. If the end is met, start from the beginning.
// Take the child and starting from the gene after the second crossover point, copy the female parent gene that is considered in the current iteration, only if it is not present yet in the child’s chromosome. If the end is met, start from the beginning.
// The operation is finished if all empty spaces in the child chromosome are filled.
// Optionally, swap the roles of female and male parents, and repeat the whole process to produce a second offspring.
// todo: Make this lower level i.e array or slice types baybeeee
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

// takes in a population and produces a children from the input
fn reproduction(population: &Vec<Solution>) -> Vec<Solution> {
    // pop need's to be even!
    assert!(population.len() % 2 == 0);
    let children = population
        .chunks(2)
        .into_iter()
        .filter_map(|chunk| {
            chunk
                .get(0)
                .and_then(|first| chunk.get(1).map(|second| (first, second)))
        })
        .map(|(a, b)| order_crossover(&a.route, &b.route))
        .map(|child| Solution { route: child })
        .collect::<Vec<Solution>>();
    children
}

// swaps: If A-B-C-D then a mutation occurs causing A-B-D-C
// should probably pass in rand FOR TESTING
fn mutate_swap<R: Rng>(route: &Vec<Location>, rnd: &mut R) -> Vec<Location> {
    let mut new_route = route.clone();
    let a = rnd.gen_range(0..route.len());
    let b = rnd.gen_range(0..route.len());
    new_route.swap(a, b);
    new_route
}
// choices:
// new-route: If we have a section of a route, i.e A -f> B find A -g-> B
// swaps: If A-B-C-D then a mutation occurs causing A-B-D-C
// for mvp i'll use swaps coz it easy bruv
fn mutation<R: Rng>(population: &Vec<Solution>, rnd: &mut R) -> Vec<Solution> {
    population
        .iter()
        .map(|solution| Solution {
            route: mutate_swap(&solution.route, rnd),
        })
        .collect()
}

fn sample_by_percentage<R: Rng>(
    population: &Vec<Solution>,
    rng: &mut R,
    percent: PercentageInteger,
) -> (Vec<Solution>, Vec<Solution>) {
    let mut sampled = population.clone();
    sampled.shuffle(rng);
    let amount_to_take = percent.apply_to(population.len());
    let (sample, remaining) = sampled.split_at(amount_to_take);
    (
        sample.into_iter().map(|s| s.to_owned()).collect(),
        remaining.into_iter().map(|s| s.to_owned()).collect(),
    )
}

// this should expose the score and a useful struct
fn next_generation(population: &Vec<Solution>) -> Vec<Solution> {
    // todo: add a seed, or probably pass the same RNG from the TL
    let mut rand = rand::thread_rng();
    let (mut current_generation, _expired) = selection(&population);
    let children = reproduction(&population);
    current_generation.extend(children);
    let (to_mutate, mut dont_mutate) =
        sample_by_percentage(&population, &mut rand, Percentage::from(10));
    dont_mutate.extend(mutation(&to_mutate, &mut rand));
    dont_mutate
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
