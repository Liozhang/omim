import re
import json
import time
from collections import defaultdict

from omim.core import OMIM


# ---------------------------------------------------------------------------
# text section IDs (in display order as they appear on OMIM pages)
# ---------------------------------------------------------------------------
TEXT_SECTION_MAP = [
    # phenotype sections
    ('Text', 'text'),
    ('Description', 'description'),
    ('ClinicalFeatures', 'clinical_features'),
    ('Inheritance', 'inheritance'),
    ('Cytogenetics', 'cytogenetics'),
    ('Mapping', 'mapping'),
    ('MolecularGenetics', 'molecular_genetics'),
    ('GenotypePhenotypeCorrelations', 'genotype_phenotype_correlations'),
    ('Pathogenesis', 'pathogenesis'),
    ('Diagnosis', 'diagnosis'),
    ('ClinicalManagement', 'clinical_management'),
    ('PopulationGenetics', 'population_genetics'),
    ('Evolution', 'evolution'),
    ('AnimalModel', 'animal_model'),
    ('History', 'history'),
    ('SeeAlso', 'see_also'),
    ('NewbornScreening', 'newborn_screening'),
    ('Nomenclature', 'nomenclature'),
    # gene sections
    ('Cloning', 'cloning'),
    ('GeneStructure', 'gene_structure'),
    ('GeneFunction', 'gene_function'),
    ('BiochemicalFeatures', 'biochemical_features'),
]


class Entry(OMIM):
    """Entry Parser For Given MIM - v2.0 with full content extraction."""

    def __init__(self, **kwarg):
        super(Entry, self).__init__(**kwarg)

    # ------------------------------------------------------------------
    # main entry point
    # ------------------------------------------------------------------

    def parse(self, mim):
        """Parse a single MIM entry page. Returns dict with all extracted data."""
        data = defaultdict(list)
        data['mim_number'] = mim

        url = self.omim_url + f'/entry/{mim}'

        # fetch page with retry
        while True:
            try:
                soup = self.get_soup(url)
                prefix = soup.select_one('#title').find_next_sibling('div').select_one('.h3 strong')
                break
            except AttributeError:
                time.sleep(3)
                self.logger.warning(f'Retrying: {url}')

        # core metadata
        data['prefix'] = prefix.text.strip() if prefix else ''
        data['title'] = self._parse_title(soup)
        data['references'] = self._parse_references(soup)

        # geneMap / phenotypeMap tables
        self._parse_gene_phenotype_maps(soup, data)

        # === NEW v2.0 ===
        data['text_sections'] = self._parse_text_sections(soup)
        data['allelic_variants'] = self._parse_allelic_variants(soup)
        data['clinical_synopsis'] = self._parse_clinical_synopsis(soup)
        data['phenotypic_series'] = self._parse_phenotypic_series(soup)

        return dict(data)

    # ------------------------------------------------------------------
    # sub-parsers
    # ------------------------------------------------------------------

    def _parse_title(self, soup):
        el = soup.select_one('#preferredTitle')
        if el:
            h3 = el.find_next_sibling('h3')
            if h3:
                return h3.text.strip()
        return ''

    def _parse_references(self, soup):
        ref = soup.select_one('#mimReferencesFold')
        if not ref:
            ref = soup.select_one('#referencesFold')
        if not ref:
            return ''
        references = re.findall(r'PubMed: (\d+)', ref.text)
        return ', '.join(references)

    def _parse_gene_phenotype_maps(self, soup, data):
        for xmap in ('phenotypeMap', 'geneMap'):
            res = soup.select_one(f'#{xmap}')
            if not res:
                continue
            table = res.parent.select_one('table')
            if not table:
                continue
            thead = table.select_one('thead')
            tbody = table.select_one('tbody')
            if not thead or not tbody:
                continue
            keys = [th.text.strip() for th in thead.select('th')]
            keys = [' '.join(key.split()).replace(' Clinical Synopses', '') for key in keys]
            for tr in tbody.select('tr'):
                row = [td.text.strip() for td in tr.select('td')]
                if len(row) == len(keys):
                    values = row
                else:
                    # handle rowspan merging
                    values = [data[xmap][-1].get(keys[0], '')] + row if data[xmap] else row
                data[xmap].append(dict(zip(keys, values)))

    # ------------------------------------------------------------------
    # v2.0: text sections
    # ------------------------------------------------------------------

    def _parse_text_sections(self, soup):
        """Extract all text subsections (Description, ClinicalFeatures, etc.)
        from the MIM page. Returns JSON string mapping section_key -> text."""
        sections = {}
        for sec_id, field_name in TEXT_SECTION_MAP:
            fold = soup.select_one(f'#mim{sec_id}Fold')
            if not fold:
                continue
            # get text, removing excessive whitespace
            text = fold.get_text(separator=' ', strip=True)
            text = re.sub(r'\s{2,}', ' ', text)
            if text:
                sections[field_name] = text
        return json.dumps(sections, ensure_ascii=False) if sections else None

    # ------------------------------------------------------------------
    # v2.0: allelic variants
    # ------------------------------------------------------------------

    def _parse_allelic_variants(self, soup):
        """Parse the allelic variants section of a gene entry.
        Returns list of variant dicts, or None if no variants section."""
        av_fold = soup.select_one('#mimAllelicVariantsFold')
        if not av_fold:
            return None

        # Each variant is anchored by <a class="mim-anchor" id="0001"></a>
        anchors = av_fold.find_all('a', class_='mim-anchor')
        if not anchors:
            return None

        # Find the container div that holds all variant blocks as children.
        # Walk up from the first anchor to find a div whose children contain
        # all anchors, then walk its immediate children grouping by anchor.
        container = anchors[0]
        while (container.parent is not None
               and container.parent.name == 'div'
               and not container.parent.get('id', '').startswith('mim')):
            container = container.parent

        # container is now the innermost div enclosing at least the first anchor.
        # Walk up one more level if needed to reach the div that holds all
        # variant blocks as direct children, not just one.
        if container.parent and container.parent.name == 'div':
            parent_anchors = container.parent.find_all('a', class_='mim-anchor')
            if len(parent_anchors) >= len(anchors):
                container = container.parent

        # Group child elements by anchor: each anchor starts a new variant block
        variants = []
        current_variant_id = None
        current_elements = []

        for child in container.find_all(recursive=False):
            anchor_in_child = (child.find('a', class_='mim-anchor')
                               if hasattr(child, 'find') else None)
            if anchor_in_child:
                vid = anchor_in_child.get('id', '').strip()
                if vid.startswith('0'):
                    # Save previous variant
                    if current_variant_id and current_elements:
                        variant = self._parse_variant_elements(
                            current_variant_id, current_elements)
                        if variant:
                            variants.append(variant)
                    current_variant_id = vid
                    current_elements = [child]
                else:
                    # Non-numeric anchor - ignore this block
                    if current_variant_id and current_elements:
                        variant = self._parse_variant_elements(
                            current_variant_id, current_elements)
                        if variant:
                            variants.append(variant)
                    current_variant_id = None
                    current_elements = []
            elif current_variant_id:
                current_elements.append(child)

        # Don't forget the last variant
        if current_variant_id and current_elements:
            variant = self._parse_variant_elements(
                current_variant_id, current_elements)
            if variant:
                variants.append(variant)

        return variants if variants else None

    def _parse_variant_elements(self, variant_id, elements):
        """Parse a variant from a list of sibling elements within the container."""
        variant = {'variant_id': '.' + variant_id}

        text_parts = []
        pubmed_ids = []

        for el in elements:
            if not hasattr(el, 'get_text'):
                continue
            text = el.get_text(separator=' ', strip=True)
            text = re.sub(r'\s{2,}', ' ', text)
            if text:
                text_parts.append(text)

            # --- phenotype name from <h4><strong>.NNNN NAME</strong></h4>
            h4 = el.find('h4') or (el if el.name == 'h4' else None)
            if h4:
                strong = h4.find('strong')
                if strong:
                    name = strong.get_text(strip=True)
                    name = re.sub(r'^\.\d+\s+', '', name)
                    name = name.rstrip(',; ')
                    variant.setdefault('phenotype_name', name)

            # --- PubMed IDs from mim-tip-reference links
            for ref in el.find_all('a', class_='mim-tip-reference'):
                pmid = ref.get('pmid', '')
                if pmid:
                    pubmed_ids.append(pmid)

        full_text = ' '.join(text_parts)
        full_text = full_text.replace('\xa0', ' ')
        variant['description'] = full_text[:10000] if full_text else None
        variant['pubmed_ids'] = ','.join(pubmed_ids) if pubmed_ids else None

        # --- extract gene_symbol, mutation, rsid, RCV from text
        for text in text_parts:
            # CFTR, PHE508DEL (rs113993960) pattern
            m = re.search(r'(\b[A-Z][A-Z0-9]+)\s*,\s*([A-Z].+?)\s*(?:\(|$)', text)
            if m and not variant.get('gene_symbol'):
                variant['gene_symbol'] = m.group(1)
                variant['mutation'] = m.group(2).strip()

            rs_match = re.search(r'(rs\d+)', text)
            if rs_match and not variant.get('rsid'):
                variant['rsid'] = rs_match.group(1)

            rcv_match = re.findall(r'(RCV\d+)', text)
            if rcv_match and not variant.get('clinvar_rcvs'):
                variant['clinvar_rcvs'] = ','.join(rcv_match[:50])

        # Ensure all keys exist
        for key in ('phenotype_name', 'gene_symbol', 'mutation', 'rsid',
                     'clinvar_rcvs', 'description', 'pubmed_ids'):
            variant.setdefault(key, None)

        return variant

    # ------------------------------------------------------------------
    # v2.0: clinical synopsis
    # ------------------------------------------------------------------

    def _parse_clinical_synopsis(self, soup):
        """Parse the clinical synopsis section (# prefix entries).
        Returns JSON string with structured phenotype data including ontology IDs."""
        cs_fold = soup.select_one('#mimClinicalSynopsisFold')
        if not cs_fold:
            return None

        # Clinical synopsis is a hierarchy of divs:
        #   .h5 strong  -> major category (e.g. INHERITANCE, CARDIOVASCULAR)
        #   .h5 em      -> subcategory   (e.g. Heart, Vascular)
        #   plain text   -> feature item with ontology IDs in .mim-feature-ids
        result = {}
        current_category = None
        current_subcategory = None

        # get top-level children divs
        container = cs_fold.find('div') or cs_fold
        for child in container.find_all('div', recursive=False):
            self._traverse_synopsis(child, result, None, None)

        return json.dumps(result, ensure_ascii=False) if result else None

    def _traverse_synopsis(self, element, result, category, subcategory):
        """Recursively traverse clinical synopsis divs."""
        # check for category heading (.h5 strong)
        strong = element.find('strong')
        em = element.find('em')

        if strong and not em:
            # major category
            category = strong.get_text(strip=True)
            if category == 'Close':
                return
            if category not in result:
                result[category] = {}
            for child in element.find_all('div', recursive=False):
                self._traverse_synopsis(child, result, category, None)
            return

        if em:
            # subcategory
            subcategory = em.get_text(strip=True)
            if category and subcategory:
                if subcategory not in result.get(category, {}):
                    result[category][subcategory] = {'items': [], 'xrefs': {}}
            for child in element.find_all('div', recursive=False):
                self._traverse_synopsis(child, result, category, subcategory)
            return

        # feature items — text with ontology IDs
        text = element.get_text(separator=' ', strip=True)
        text = re.sub(r'\s{2,}', ' ', text)
        text = text.strip('- ')

        if not text:
            return

        # extract ontology IDs from .mim-feature-ids spans
        xrefs = {}
        for feat_span in element.find_all('span', class_='mim-feature-ids'):
            for link in feat_span.find_all('a'):
                href = link.get('href', '')
                txt = link.get_text(strip=True)
                # determine ontology source from href
                if 'HP:0' in txt:
                    xrefs.setdefault('HPO', []).append(txt.split(':')[-1] if ':' in txt else txt)
                elif 'SNOMEDCT' in href or 'SNOMEDCT' in feat_span.get_text():
                    xrefs.setdefault('SNOMEDCT', []).append(txt)
                elif 'ICD10CM' in href:
                    xrefs.setdefault('ICD10CM', []).append(txt)
                elif 'ICD9CM' in href:
                    xrefs.setdefault('ICD9CM', []).append(txt)
                elif 'UMLS:' in feat_span.get_text():
                    xrefs.setdefault('UMLS', []).append(txt)

        # store the feature item
        if category and subcategory:
            bucket = result[category].setdefault(subcategory, {'items': [], 'xrefs': {}})
            if text not in bucket['items']:
                bucket['items'].append(text)
            for src, ids in xrefs.items():
                for id_ in ids:
                    if id_ not in bucket['xrefs'].get(src, []):
                        bucket['xrefs'].setdefault(src, []).append(id_)
        elif category:
            bucket = result[category].setdefault('_general', {'items': [], 'xrefs': {}})
            if text not in bucket['items']:
                bucket['items'].append(text)
            for src, ids in xrefs.items():
                for id_ in ids:
                    if id_ not in bucket['xrefs'].get(src, []):
                        bucket['xrefs'].setdefault(src, []).append(id_)

    # ------------------------------------------------------------------
    # v2.0: phenotypic series
    # ------------------------------------------------------------------

    def _parse_phenotypic_series(self, soup):
        """Extract phenotypic series MIM numbers.
        Returns comma-separated string or None."""
        ps_toggle = soup.select_one('#mimPhenotypicSeriesToggle')
        if not ps_toggle:
            return None

        # Go up to the containing div and find MIM number links
        container = ps_toggle.parent
        while container and container.name != 'div':
            container = container.parent

        if not container:
            return None

        mim_list = []
        for link in container.find_all('a'):
            href = link.get('href', '')
            # Look for /entry/NNNNNN links
            m = re.search(r'/entry/(\d+)', href)
            if m:
                mim = m.group(1)
                if mim not in mim_list:
                    mim_list.append(mim)

        return ','.join(mim_list) if mim_list else None


if __name__ == '__main__':
    from pprint import pprint

    entry = Entry()

    # * gene entries
    # data = entry.parse('612367')    # one geneMap
    # data = entry.parse('607093')    # MTHFR - geneMap + allelicVariants
    data = entry.parse('602421')    # CFTR - many allelicVariants
    # data = entry.parse('300050')    # no geneMap
    # data = entry.parse('109690')    # multiple geneMap

    # # phenotype entries
    # data = entry.parse('219700')    # CF - clinicalSynopsis + phenotypicSeries

    # other types
    # data = entry.parse('100500')    # moved
    # data = entry.parse('618428')    # removed
    # data = entry.parse('100650')    # +
    # data = entry.parse('100070')    # %
    # data = entry.parse('100100')    # #
    # data = entry.parse('100050')    # other

    pprint(data)
