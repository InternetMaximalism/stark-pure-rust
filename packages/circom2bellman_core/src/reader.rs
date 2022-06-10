use crate::r1csfile::*;
use bytes::Buf;

pub fn read_bytes(_bytes: &[u8]) -> R1csContents {
    let mut p = _bytes;
    let magic = p.get_u32_le();
    assert_eq!(u32::from_le_bytes(*b"r1cs"), magic);
    let version = Version(p.get_u32_le());
    assert_eq!(version, Version(1));
    let n_section = NSection(p.get_u32_le());
    assert_eq!(3, n_section.0);

    let section_type_header = p.get_u32_le();
    assert_eq!(SectionType::HeaderSection as u32, section_type_header);
    let _section_size = p.get_u64_le();
    let field_size = p.get_u32_le();
    let mut prime_number = [0u8; 32]; // little endian
    for i in 0..32 {
        prime_number[i] = p.get_u8();
    }
    let n_wires = p.get_u32_le();
    let n_public_outputs = p.get_u32_le();
    let n_public_inputs = p.get_u32_le();
    let n_private_inputs = p.get_u32_le();
    let n_labels = p.get_u64_le();
    let n_constraints = p.get_u32_le();
    let header = Header {
        field_size,
        prime_number,
        n_wires,
        n_public_outputs,
        n_public_inputs,
        n_private_inputs,
        n_labels,
        n_constraints,
    };

    let section_type_constraints = p.get_u32_le();
    assert_eq!(
        SectionType::ConstraintSection as u32,
        section_type_constraints
    );
    let _section_size = p.get_u64_le();
    let mut constraint_vectors = Vec::<Constraint>::new();
    for _ in 0..(n_constraints as usize) {
        let mut factors: Vec<Factor> = Vec::with_capacity(3);
        for _ in 0..3 {
            let n_coefficient = p.get_u32_le();
            let mut coefficients = Vec::<Coefficient>::new();
            for _ in 0..n_coefficient {
                let wire_id = p.get_u32_le();
                let mut value = [0u8; 32]; // little endian
                for i in 0..32 {
                    value[i] = p.get_u8();
                }
                let coefficient = Coefficient { wire_id, value };
                coefficients.push(coefficient)
            }
            let factor = Factor {
                n_coefficient,
                coefficients,
            };
            factors.push(factor);
        }
        let constraint = Constraint { factors };
        constraint_vectors.push(constraint)
    }
    let constraints = Constraints(constraint_vectors);

    // let section_type_map = p.get_u32_le();
    // assert_eq!(SectionType::Wire2LabelIdMapSection, section_type_map);
    // p.get_u64_le();
    // let mut label_ids = Vec::<LabelId>::new();
    // println!("{}", n_labels);
    // for i in 0..(n_labels as usize) {
    //     println!("10000");
    //     let label_id = LabelId(p.get_u64_le());
    //     println!("20000");
    //     label_ids.push(label_id);
    // }
    // let wire2_label_id_map = Wire2LabelledMap(label_ids);

    
    R1csContents {
        version,
        header,
        constraints,
        // wire2_label_id_map
    }
}
