use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_pwmantel::{Alternative, LabeledMatrix, Method, pwmantel, write_results};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-pwmantel", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Two or more distance matrices (lsmat TSV: blank corner, id header, rows).
    #[arg(required = true, num_args = 2..)]
    matrices: Vec<PathBuf>,

    #[arg(short = 'm', long, default_value = "pearson")]
    method: String,

    #[arg(short = 'p', long, default_value_t = 999)]
    permutations: usize,

    #[arg(short = 'a', long, default_value = "two-sided")]
    alternative: String,

    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let method = Method::parse(&self.method).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "invalid method '{}' (pearson|spearman)",
                self.method
            ))
        })?;
        let alternative = Alternative::parse(&self.alternative).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "invalid alternative '{}' (two-sided|greater|less)",
                self.alternative
            ))
        })?;
        let seed = self.common.seed_rng();

        let dms = self
            .matrices
            .iter()
            .map(|p| LabeledMatrix::read(open(p)?, label_of(p)))
            .collect::<Result<Vec<_>>>()?;

        let results = pwmantel(&dms, method, self.permutations, alternative, seed)?;

        let mut out: Box<dyn Write> = if self.output == "-" && self.common.json {
            Box::new(BufWriter::new(std::io::sink()))
        } else if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        write_results(&mut out, &results, method, self.permutations, alternative)?;
        out.flush().map_err(RsomicsError::Io)?;

        if !self.common.quiet {
            eprintln!(
                "pwmantel: {} method, {} pairs over {} matrices",
                method.name(),
                results.len(),
                dms.len()
            );
        }
        Ok(())
    }
}

fn label_of(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

fn open(path: &Path) -> Result<BufReader<File>> {
    File::open(path)
        .map(BufReader::new)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Pairwise Mantel test across N distance matrices.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.stats.distance.pwmantel",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("PMID:6018555"),
    }),
    usage_lines: &[
        "<dm1.tsv> <dm2.tsv> [dm3.tsv ...] [-m pearson] [-p 999] [-a two-sided] [-o out.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('m'),
                long: "method",
                aliases: &[],
                value: Some("<pearson|spearman>"),
                type_hint: Some("String"),
                required: false,
                default: Some("pearson"),
                description: "Correlation method; spearman ranks the distances first.",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "permutations",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("999"),
                description: "Permutations for the p-value; 0 skips it (p = nan).",
                why_default: None,
            },
            FlagSpec {
                short: Some('a'),
                long: "alternative",
                aliases: &[],
                value: Some("<two-sided|greater|less>"),
                type_hint: Some("String"),
                required: false,
                default: Some("two-sided"),
                description: "Alternative hypothesis for the permutation p-value.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Pairwise Pearson Mantel over three matrices",
            command: "rsomics-pwmantel a.tsv b.tsv c.tsv -o results.tsv",
        },
        Example {
            description: "Spearman, fixed seed",
            command: "rsomics-pwmantel a.tsv b.tsv c.tsv -m spearman --seed 42",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
