#[derive(Debug, PartialEq)]
pub struct R1csContents {
    pub version: Version,
    pub header: Header,
    pub constraints: Constraints,
    // pub wire2_label_id_map: Wire2LabelledMap
}

#[derive(Debug, PartialEq)]
pub struct Version(pub u32);

#[derive(Debug, PartialEq)]
pub struct NSection(pub u32);

#[derive(Debug, PartialEq)]
pub enum SectionType {
    Other,
    HeaderSection,
    ConstraintSection,
    // Wire2LabelIdMapSection,
}

#[derive(Debug, PartialEq)]
pub struct SectionSize(u64);

#[derive(Debug, PartialEq)]
pub struct Header {
    pub field_size: u32,
    pub prime_number: [u8; 32],
    pub n_wires: u32,
    pub n_public_outputs: u32,
    pub n_public_inputs: u32,
    pub n_private_inputs: u32,
    pub n_labels: u64,
    pub n_constraints: u32,
}

#[derive(Debug, PartialEq)]
pub struct Constraints(pub Vec<Constraint>);

#[derive(Debug, PartialEq)]
pub struct Constraint {
    pub factors: Vec<Factor>,
}

#[derive(Debug, PartialEq)]
pub struct Factor {
    pub n_coefficient: u32,
    pub coefficients: Vec<Coefficient>,
}

#[derive(Debug, PartialEq)]
pub struct Coefficient {
    pub wire_id: u32,
    pub value: [u8; 32],
}

/*
#[derive(Debug,PartialEq)]
pub struct Wire2LabelledMap(pub Vec<LabelId>);
*/

#[derive(Debug, PartialEq)]
pub struct LabelId(pub u64);

pub trait VerifyForm {
    fn verify_form(&self) -> bool;
}

impl VerifyForm for R1csContents {
    fn verify_form(&self) -> bool {
        let is_valid_constraints = self.header.n_constraints as usize == self.constraints.0.len();
        // let is_valid_labels = self.header.n_labels as usize == self.wire2_label_id_map.0.len();
        let mut is_valid_factors = true;
        for constraint in self.constraints.0.iter() {
            if constraint.factors.len() != 3
                || !constraint.factors[0].verify_form()
                || !constraint.factors[1].verify_form()
                || !constraint.factors[2].verify_form()
            {
                is_valid_factors = false;
                break;
            }
        }
        return is_valid_constraints && /*is_valid_labels &&*/ is_valid_factors;
    }
}

impl VerifyForm for Factor {
    fn verify_form(&self) -> bool {
        self.n_coefficient as usize == self.coefficients.len()
    }
}
