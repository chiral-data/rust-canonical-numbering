use crate::core;
use crate::core::graph::*;
use super::atom;
use super::extendable_hash;
use super::local_symmetry;
use super::molecule;

pub fn reduced_edges(
    atoms: &Vec<atom::Atom>, 
    atom_indexes: &Vec<usize>
) -> Vec<(usize, usize, usize)> {
    atom_indexes.iter()
        .map(|&ai| atoms[ai].bonds.iter().map(move |b| (ai, b.tid, b.fixed_hash_value())))
        .flatten()
        .collect()
}

pub fn canonical_numbering_and_symmetry_perception(
    atoms: &Vec<atom::Atom>,
    orbits_after_partition: &mut Vec<core::orbit_ops::Orbit>,
    orbits_symmetry: &mut Vec<core::orbit_ops::Orbit>,
    numbering: &mut Vec<usize>,
) {
    let indexes_all: Vec<usize> = (0..atoms.len()).collect();
    let vv = core::graph::VertexVec::init(indexes_all, atoms.to_vec()); 
    numbering.clear();
    // calculate and save the givp result for comparison
    let givp_start = std::time::Instant::now();
    core::givp::run::<extendable_hash::AtomExtendable>(&vv, numbering, orbits_after_partition);
    let givp_duration = givp_start.elapsed();
    let count_of_atoms = vv.len();
    if cfg!(debug_assertions) {
        // println!("numbering: {:?}", numbering);
    }

    let mut rg = core::reduce::reducible_graph::ReducibleGraph {
        vv: vv,
        mapping: vec![],
        boundary_edges: vec![],
        orbits_after_partition: orbits_after_partition.clone(),
        numbering: numbering.clone() 
    };

    let cnap_start = std::time::Instant::now();
    core::symmetry_perception::<extendable_hash::AtomExtendable>(&mut rg, orbits_symmetry, reduced_edges, local_symmetry::get_local_symmetric_orbits, 200);
    let cnap_duration = cnap_start.elapsed();
    if cfg!(debug_assertions) {
        let cnap_required: bool = orbits_after_partition.len() != 0;
        println!("atom_count\tgivp_duration\tcnap_required\tcnap_duration");
        println!("{:?}\t{:?}\t{:?}\t{:?}", count_of_atoms, givp_duration, cnap_required, cnap_duration);
    }
}

fn symmetry_perception_givp(
    smiles: &String
) -> (Vec<core::orbit_ops::Orbit>, Vec<usize>, core::graph::VertexVec<atom::Atom>) {
    let mol = molecule::Molecule::from_smiles(&smiles);
    if cfg!(debug_assertions) {
        println!("{}", mol.smiles_with_index(&smiles, &vec![]));
    }
    let indexes_all: Vec<usize> = (0..mol.atoms.len()).collect();
    let vv = core::graph::VertexVec::init(indexes_all, mol.atoms); 
    let mut numbering: Vec<usize> = vec![];
    let mut orbits_after_partition: Vec<core::orbit_ops::Orbit> = vec![];
    core::givp::run::<extendable_hash::AtomExtendable>(&vv, &mut numbering, &mut orbits_after_partition);

    (orbits_after_partition, numbering, vv)
}

fn symmetry_perception_cnap(
    smiles: &String
) -> Vec<core::orbit_ops::Orbit> {
    let (orbits_givp, numbering, vv) = symmetry_perception_givp(&smiles);
    let mut orbits_cnap: Vec<core::orbit_ops::Orbit> = vec![];

    if orbits_givp.len() != 0 {
        let mut rg = core::reduce::reducible_graph::ReducibleGraph {
            vv: vv,
            mapping: vec![],
            boundary_edges: vec![],
            orbits_after_partition: orbits_givp.clone(),
            numbering: numbering
        };
        core::symmetry_perception::<extendable_hash::AtomExtendable>(&mut rg, &mut orbits_cnap, reduced_edges, local_symmetry::get_local_symmetric_orbits, 200);
    }

    core::orbit_ops::orbits_sort(&mut orbits_cnap);
    orbits_cnap
}

fn canonical_numbering(
    smiles: &String,
) -> Vec<usize> {
    let mol = molecule::Molecule::from_smiles(&smiles);
    if cfg!(debug_assertions) {
        println!("{}", mol.smiles_with_index(&smiles, &vec![]));
    }
    let indexes_all: Vec<usize> = (0..mol.atoms.len()).collect();
    let mut vv = core::graph::VertexVec::init(indexes_all, mol.atoms.clone()); 
    let mut numbering: Vec<usize> = vec![];
    let mut orbits_givp: Vec<core::orbit_ops::Orbit> = vec![];
    core::givp::run::<extendable_hash::AtomExtendable>(&vv, &mut numbering, &mut orbits_givp);

    println!("before break-tie: \n{}", mol.smiles_with_index(&smiles, &numbering));
    println!("numbering: {:?}, {:?}", numbering, orbits_givp);

    while orbits_givp.len() > 0 {
        let mut symmetrical_atoms: Vec<usize> = orbits_givp.iter()
            .map(|orbit| orbit[0])
            .collect();
        symmetrical_atoms.sort_by_key(|&idx| numbering[idx]);
        symmetrical_atoms.reverse();
        if let Some(vi) = vv.all_vertices_mut().get_mut(symmetrical_atoms[0]) {
            vi.break_symmetry_vertex();
        }
        let fixed_hash_values: Vec<core::graph::VertexFixedHashValue> = vv.all_vertices().iter()
            .map(|v| v.fixed_hash_value())
            .collect();
        core::givp::partition_vertices_with::<extendable_hash::AtomExtendable>(&vv.valid_indexes(), vv.all_vertices(), &fixed_hash_values, &mut numbering, &mut orbits_givp);
    }

    println!("after break-tie: \n{}", mol.smiles_with_index(&smiles, &numbering));

    numbering
}

#[cfg(test)]
mod test_molecule_workflow {
    use super::*;
    use crate::ext::molecule;

    #[test]
    fn test_examples() {
        type InputType1 = String;
        let test_data: Vec<InputType1> = vec![
            // 
            // *** SOLVED ***
            // "CCn1c2ccc3cc2c2cc(ccc21)C(=O)c1ccc(cc1)Cn1cc[n+](c1)Cc1ccc(cc1)-c1cccc(c1C(=O)O)-c1ccc(cc1)C[n+]1ccn(c1)Cc1ccc(cc1)C3=O", // chembl 15,
            // "CC(C)(CCCOc1cc(Cl)c(OCCCC(C)(C)C(=O)O)cc1Cl)C(=O)O", // 4631
            // "C[N+](C)(CCCCCC[N+](C)(C)CCCN1C(=O)C2C3c4ccccc4C(c4ccccc43)C2C1=O)CCCN1C(=O)c2ccccc2C1=O", // 6053 separable graph 
            // "N[C@@H](Cc1cnc(C23CC4CC(CC(C4)C2)C3)[nH]1)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1cnc(C23CC4CC(CC(C4)C2)C3)[nH]1)C(=O)NCc1ccccc1", // 7844 separable graph 
            // "OCCCCCNCc1c2ccccc2c(CNCCCCCO)c2ccccc12", // 23218 
            // // "NC[C@@H]1O[C@H](O[C@@H]2[C@@H](CSCCNC(=S)NCCCCn3c(=O)c4ccc5c6ccc7c(=O)n(CCCCNC(=S)NCCSC[C@H]8O[C@@H](O[C@@H]9[C@@H](O)[C@H](N)C[C@H](N)[C@H]9O[C@H]9O[C@H](CN)[C@@H](O)[C@H](O)[C@H]9N)[C@H](O)[C@@H]8O[C@H]8O[C@@H](CN)[C@@H](O)[C@H](O)[C@H]8N)c(=O)c8ccc(c9ccc(c3=O)c4c59)c6c78)O[C@@H](O[C@@H]3[C@@H](O)[C@H](N)C[C@H](N)[C@H]3O[C@H]3O[C@H](CN)[C@@H](O)[C@H](O)[C@H]3N)[C@@H]2O)[C@H](N)[C@@H](O)[C@@H]1O", // 52881 
            // "CC1(C)c2ccc([nH]2)C2(C)CCCCNC(=O)c3cccc(n3)C(=O)NCCCCC(C)(c3ccc1[nH]3)c1ccc([nH]1)C(C)(C)c1ccc2[nH]1", // 4971 interesting example, 8 vertices cycle, 2 folded symmetry
            // "O=C1NNC(=O)c2ccccc2SSc2ccccc2C(=O)NNC(=O)c2ccccc2SSc2ccccc21", // 140635
            // "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 168272
            // "O=P1([O-])OC2C3OP(=O)([O-])OP(=O)([O-])OC3C3OP(=O)([O-])OP(=O)([O-])OC3C2OP(=O)([O-])O1", // 171007
            // "C1CC1N1CN2c3nonc3N3CN(C4CC4)CN4c5nonc5N(C1)C2C34", // 199821
            // "O=P1(O)OC2C3OP(=O)(O)OP(=O)(O)OC3C3OP(=O)(O)OP(=O)(O)OC3C2OP(=O)(O)O1", // 208361
            // // "CC[n+]1ccc(-c2cc[n+](Cc3cc(C[n+]4ccc(-c5cc[n+](CC)cc5)cc4)cc(C[n+]4ccc(-c5cc[n+](Cc6cc(C[n+]7ccc(-c8cc[n+](Cc9cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)c9)cc8)cc7)cc(-[n+]7ccc(-c8cc[n+](-c9cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)c9)cc8)cc7)c6)cc5)cc4)c3)cc2)cc1", // 826428 long givp time
            // // "CC[n+]1ccc(-c2cc[n+](Cc3cc(C[n+]4ccc(-c5cc[n+](CC)cc5)cc4)cc(C[n+]4ccc(-c5cc[n+](Cc6cc(C[n+]7ccc(-c8cc[n+](Cc9cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](CC)cc%11)cc%10)c9)cc8)cc7)cc(-[n+]7ccc(-c8cc[n+](-c9cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)cc(C[n+]%10ccc(-c%11cc[n+](Cc%12cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)cc(C[n+]%13ccc(-c%14cc[n+](CC)cc%14)cc%13)c%12)cc%11)cc%10)c9)cc8)cc7)c6)cc5)cc4)c3)cc2)cc1", // 1246825
            // "BrC1CCC(Br)C(Br)CCC(Br)C(Br)CCC1Br", // 377203
            // "C[N+]1(C)CC23c4c5c6c7c8c4c4c2c2c9c%10c%11c%12c%13c9c9c%14c%15c%16c%17c%18c%19c(c8c%17c4c%16c29)C7C2c4c-%19c7c8c9c(c%14c%13c%13c9c9c8c4c4c2c6c2c5c(c=%11c5c2c4c9c5c%12%13)C%103C1)C%15C%187", // CHEMBL415840 failure case in Schneider paper
            // "CCC[C@H]1CC[C@H]([C@H]2CC[C@H](OC(=O)[C@H]3[C@@H](c4ccc(O)cc4)[C@H](C(=O)O[C@H]4CC[C@H]([C@H]5CC[C@H](CCC)CC5)CC4)[C@@H]3c3ccc(O)cc3)CC2)CC1", // CHEMBL2348759, failure case in Schneider paper
            // 
            // *** NOT SOLVED ***
            // "OC(c1ccccc1)C1(c2ccccc2)C23c4c5c6c7c8c9c(c%10c%11c2c2c4c4c%12c5c5c6c6c8c8c%13c9c9c%10c%10c%11c%11c2c2c4c4c%12c%12c5c5c6c8c6c8c%13c9c9c%10c%10c%11c2c2c4c4c%12c5c6c5c8c9c%10c2c45)C731", // 408840 beneze ball
            // "O=C(CCCc1ccc(C2(c3ccccc3)C34c5c6c7c8c9c%10c(c%11c%12c3c3c5c5c%13c6c6c7c7c9c9c%14c%10c%10c%11c%11c%12c%12c3c3c5c5c%13c%13c6c6c7c9c7c9c%14c%10c%10c%11c%11c%12c3c3c5c5c%13c6c7c6c9c%10c%11c3c56)C824)cc1)NC(CO)(CO)CO", // 267348 beneze ball
            // "O=C(CCCc1ccc(C2(c3ccccc3)C34c5c6c7c8c9c%10c(c%11c%12c3c3c5c5c%13c6c6c7c7c9c9c%14c%10c%10c%11c%11c%12c%12c3c3c5c5c%13c%13c6c6c7c9c7c9c%14c%10c%10c%11c%11c%12c3c3c5c5c%13c6c7c6c9c%10c%11c3c56)C824)cc1)NC(CO)(CO)CO", // 267348
            // r#"C[C@H](CC[C@@H]([C@@H]([C@H](C)C[C@H](C(=C)/C(=C/CO)/C)O)O)OS(=O)(=O)[O-])[C@H]([C@@H](C)[C@H]1[C@@H]([C@@H]([C@H]2[C@H](O1)[C@@H](C[C@]3([C@H](O2)C[C@H]4[C@H](O3)C[C@]5([C@H](O4)[C@H]([C@H]6[C@H](O5)C[C@H]([C@H](O6)[C@@H]([C@H](C[C@H]7[C@@H]([C@@H]([C@H]8[C@H](O7)C[C@H]9[C@H](O8)C[C@H]1[C@H](O9)[C@H]([C@@H]2[C@@H](O1)[C@@H]([C@H]([C@@H](O2)[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@@H](C[C@H](C[C@H]1[C@@H]([C@H]([C@H]2[C@@H](O1)C[C@H]([C@@H](O2)[C@H]1[C@@H](C[C@]2([C@H](O1)[C@@H]([C@]1([C@H](O2)C[C@]2([C@H](O1)CC[C@]1([C@H](O2)C[C@]2([C@H](O1)C[C@H]1[C@H](O2)CC[C@H](O1)[C@]1([C@@H](C[C@H]2[C@](O1)(C[C@H]1[C@](O2)(CC[C@]2([C@H](O1)C[C@H]1[C@](O2)(C[C@H]2[C@H](O1)C/C=C\[C@H]1[C@H](O2)C[C@H]2[C@](O1)(C[C@]1([C@H](O2)C[C@H]2[C@](O1)(CC[C@H](O2)[C@H]([C@@H](C[C@@H](C)[C@@H](C)CC=C)O)O)C)C)C)C)C)C)C)O)C)C)C)C)C)O)C)O)O)O)O)O)O)O)O)O)O)O)O)O)OS(=O)(=O)[O-])O)O)O)O)C)C)O)O)O)O"#, // Maitotoxin
            // "OC(=O)c1cc2Cc3cc(Cc4cc(Cc5cc(Cc(c2)c1)cc(c5)C(O)=O)cc(c4)C(O)=O)cc(c3)C(O)=O", // graph reduction demo
            "C1C2CC3CC1CC(C2)C3", // example from nauty, https://pallini.di.uniroma1.it/Introduction.html
       ].into_iter().map(|s| s.to_string()).collect();

        for td in test_data.iter() {
            let smiles = td.clone();
            let mol = molecule::molecule::Molecule::from_smiles(&smiles);
            if cfg!(debug_assertions) {
                println!("{}", mol.smiles_with_index(&smiles, &vec![]));
            }

            let mut orbits_partitioned: Vec<core::orbit_ops::Orbit> = vec![];
            let mut orbits_symmetry: Vec<core::orbit_ops::Orbit> = vec![];
            let mut numbering: Vec<usize> = vec![];
            molecule::workflow::canonical_numbering_and_symmetry_perception(&mol.atoms, &mut orbits_partitioned, &mut orbits_symmetry, &mut numbering);
            println!("{}", mol.smiles_with_index(&smiles, &numbering));
            if cfg!(debug_assertions) {
                core::orbit_ops::orbits_sort(&mut orbits_partitioned);
                core::orbit_ops::orbits_sort(&mut orbits_symmetry);
                println!("GIVP: {:?}\nCNAP: {:?}", orbits_partitioned, orbits_symmetry);
            }
            assert_eq!(core::orbit_ops::orbits_equal(&orbits_partitioned, &orbits_symmetry), true);
        }
    }

    #[test]
    fn test_symmetry_perception() {
        type InputType1 = String;
        type ReturnType1 = Vec<core::orbit_ops::Orbit>;
        type ReturnType2 = Vec<core::orbit_ops::Orbit>;
        type ReturnType3 = Vec<usize>;
        let test_data: Vec<(InputType1, ReturnType1, ReturnType2, ReturnType3)> = vec![
            (
                "C(C)(C)CCN",
                vec![vec![1, 2]],
                vec![vec![1, 2]],
                vec![6, 2, 1, 5, 4, 3],
            ),
            (
                "C(C)(C)CCNCCC(C)(C)",
                vec![vec![0, 8], vec![1, 2, 9, 10], vec![3, 7], vec![4, 6]],
                vec![vec![0, 8], vec![1, 2, 9, 10], vec![3, 7], vec![4, 6]],
                vec![]
            )
        ].into_iter().map(|s| (s.0.to_string(), s.1, s.2, s.3)).collect();

        for td in test_data.iter() {
            let (smiles, orbits_givp, orbits_cnap, canon_numbering) = td.clone();
            assert_eq!(symmetry_perception_givp(&smiles).0, orbits_givp);
            assert_eq!(symmetry_perception_cnap(&smiles), orbits_cnap);
            assert_eq!(canonical_numbering(&smiles), canon_numbering);
        }
    }
}