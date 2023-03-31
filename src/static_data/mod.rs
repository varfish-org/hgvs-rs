//! Static data.

use std::io::Read;

use enum_map::{enum_map, Enum, EnumMap};
use flate2::read::GzDecoder;
use serde::Deserialize;

const GRCH37_JSON_GZ: &[u8] = include_bytes!("_data/GRCh37.json.gz");
const GRCH37_P10_JSON_GZ: &[u8] = include_bytes!("_data/GRCh37.p10.json.gz");
const GRCH38_JSON_GZ: &[u8] = include_bytes!("_data/GRCh38.json.gz");

#[derive(Debug, Deserialize, Enum, Clone, Copy, PartialEq, Eq)]
pub enum Assembly {
    Grch37,
    Grch37p10,
    Grch38,
}

impl Assembly {
    /// Deserialize assembly info from embedded compressed JSON.
    fn load_assembly_info(&self) -> AssemblyInfo {
        let payload = match self {
            Assembly::Grch37 => GRCH37_JSON_GZ,
            Assembly::Grch37p10 => GRCH37_P10_JSON_GZ,
            Assembly::Grch38 => GRCH38_JSON_GZ,
        };
        let mut d = GzDecoder::new(payload);
        let mut grch37_json = String::new();
        d.read_to_string(&mut grch37_json)
            .expect("should not happen; invalid gzip in embedded data");
        serde_json::from_str::<AssemblyInfo>(&grch37_json)
            .expect("should not happen; invalid JSON in embedded data")
    }
}

#[derive(Debug, Deserialize)]
pub struct Sequence {
    pub aliases: Vec<String>,
    pub assembly_unit: String,
    pub genbank_ac: String,
    pub length: usize,
    pub name: String,
    pub refseq_ac: String,
    pub relationship: String,
    pub sequence_role: String,
}

#[derive(Debug, Deserialize)]
pub struct AssemblyInfo {
    pub date: String,
    pub description: String,
    pub genbank_ac: String,
    pub name: String,
    pub refseq_ac: String,
    pub sequences: Vec<Sequence>,
    pub submitter: String,
}

lazy_static::lazy_static! {
    /// Provide information about the assemblies.
    pub static ref ASSEMBLY_INFOS: EnumMap<Assembly, AssemblyInfo> = enum_map! {
        Assembly::Grch37 => Assembly::Grch37.load_assembly_info(),
        Assembly::Grch37p10 => Assembly::Grch37p10.load_assembly_info(),
        Assembly::Grch38 => Assembly::Grch38.load_assembly_info(),
    };
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::static_data::{Assembly, ASSEMBLY_INFOS};

    #[test]
    fn smoke() {
        assert_eq!(ASSEMBLY_INFOS[Assembly::Grch37].sequences.len(), 92);
        assert_eq!(ASSEMBLY_INFOS[Assembly::Grch37p10].sequences.len(), 275);
        assert_eq!(ASSEMBLY_INFOS[Assembly::Grch38].sequences.len(), 455);
    }
}

// <LICENSE>
// Copyright 2023 hgvs-rs Contributors
// Copyright 2014 Bioutils Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// </LICENSE>
