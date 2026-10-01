pub enum FramingMode {
    Mdl,
    Metadata,
}

impl FramingMode {
    pub fn prefix(&self) -> &'static str {
        match self {
            Self::Mdl => "M  END",
            Self::Metadata => "$$$$",
        }
    }
}
