from sqlalchemy import (Column, Integer, Float, DECIMAL, String, Text,
                        DATETIME, ForeignKey, BOOLEAN, Index, DATE)
from sqlalchemy.orm import relationship
from sqlalchemy.orm.state import InstanceState
from sqlalchemy.ext.declarative import declarative_base


# 创建对象的基类:
Base = declarative_base()


class OMIM_DATA(Base):
    __tablename__ = 'omim'

    mim_number = Column(String(10), primary_key=True, comment='MIM Number')

    prefix = Column(String(1), comment='The prefix symbol')
    title = Column(String(300), comment='The title')
    references = Column(Text, comment='The references (comma-separated PubMed IDs)')

    geneMap = Column(Text, comment='The geneMap data (JSON)')
    phenotypeMap = Column(Text, comment='The phenotypeMap data (JSON)')

    mim_type = Column(String(20), comment='The mim_type')
    entrez_gene_id = Column(String(20), comment='The entrez_gene_id')
    ensembl_gene_id = Column(String(20), comment='The ensembl_gene_id')
    hgnc_gene_symbol = Column(String(20), comment='The hgnc_gene_symbol')

    generated = Column(DATETIME, comment='The generated time')

    # === NEW COLUMNS (v2.0) ===
    text_sections = Column(Text, comment='Full text subsections (JSON: section_name -> text)')
    clinical_synopsis = Column(Text, comment='Clinical synopsis structured data (JSON with ontology IDs)')
    phenotypic_series = Column(Text, comment='Phenotypic series MIM numbers (comma-separated)')
    parser_version = Column(String(10), comment='Parser version used to scrape this entry')

    __table_args__ = (
        Index('search_by_gene', 'hgnc_gene_symbol'),
        Index('search_by_title', 'title'),
    )

    @property
    def as_dict(self):
        return {k: v for k, v in self.__dict__.items()
                if not isinstance(v, InstanceState)}

    def __str__(self):
        return '[{mim_number} - {title}]'.format(**self.__dict__)

    __repr__ = __str__


class OMIM_ALLELIC_VARIANT(Base):
    """Allelic variants parsed from individual MIM gene entry pages."""
    __tablename__ = 'omim_allelic_variants'

    id = Column(Integer, primary_key=True, autoincrement=True)
    mim_number = Column(String(10), ForeignKey('omim.mim_number'), nullable=False,
                        comment='Parent MIM entry number')
    variant_id = Column(String(10), nullable=False,
                        comment='Variant identifier (.0001, .0002, ...)')
    phenotype_name = Column(String(300), comment='Phenotype name for this variant')
    gene_symbol = Column(String(30), comment='Gene symbol')
    mutation = Column(Text, comment='Mutation description (e.g. PHE508DEL)')
    rsid = Column(String(30), comment='dbSNP rsID')
    clinvar_rcvs = Column(Text, comment='ClinVar RCV accessions (comma-separated)')
    description = Column(Text, comment='Variant description text')
    pubmed_ids = Column(Text, comment='PubMed IDs cited in variant description')

    __table_args__ = (
        Index('idx_av_mim', 'mim_number'),
        Index('idx_av_variant_id', 'mim_number', 'variant_id'),
        Index('idx_av_rsid', 'rsid'),
        Index('idx_av_symbol', 'gene_symbol'),
    )

    def __str__(self):
        return f'[{self.mim_number}.{self.variant_id} - {self.phenotype_name}]'

    __repr__ = __str__


OMIM_DATA_COLUMNS = dict(OMIM_DATA.metadata.tables['omim'].columns)
