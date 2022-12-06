import os
import urllib.request
import urllib.parse
import urllib.error
import json
import shutil
import tempfile
import argparse

import pandas as pd
import pyranges as pr
#import pyBigWig


def http_get(url, params=None, headers=None, text=False, retry_count=10, retry_interval=1):
    # set params
    if params is not None:
        url = url + '?' + urllib.parse.urlencode(params)
    if headers is None:
        if text:
            headers = {'Accept': 'text/plain'}
        else:
            headers = {'Accept': 'application/json'}
    # main
    req = urllib.request.Request(url, headers=headers, method='GET')
    return http_send_request(req, text, retry_count, retry_interval)


def http_send_request(req, text, retry_count, retry_interval):
    n_try = 0
    while True:
        n_try += 1
        try:
            with urllib.request.urlopen(req) as response:
                if text:
                    result = response.read().decode('utf-8')
                else:
                    result = json.loads(response.read())
        except urllib.error.HTTPError as exc:
            if exc.code == 500:  # internal server error
                if n_try > retry_count:
                    print(exc.read().decode('utf-8'))
                    raise Exception(f'Exceeded maximum retry count.') from exc
                else:
                    time.sleep(retry_interval)
                    continue
            else:
                print(exc.read().decode('utf-8'))
                raise
        else:
            break

    return result


def download(url, path):
    with urllib.request.urlopen(url) as response:
        with open(path, 'wb') as outfile:
            shutil.copyfileobj(response, outfile)


############################


ENCODE_PREFIX = 'https://www.encodeproject.org'
ALL_cCREs_URL = 'https://downloads.wenglab.org/V3/GRCh38-cCREs.bed'


def query_ENCODE(url_part, params=None):
    return http_get(os.path.join(ENCODE_PREFIX, url_part), params=params, text=False)


def get_annotation_search_options():
    result = dict()
    query_result = query_ENCODE('search', params={'type': 'Annotation'})
    for item in query_result['facets']:
        result[item['field']] = [x['key'] for x in item['terms']]

    return result


def get_available_organs():
    return sorted(get_annotation_search_options()['biosample_ontology.organ_slims'])


def get_all_cCREs_gr():
    fd, tmpfile_path = tempfile.mkstemp(suffix='.bed', dir=os.getcwd())
    os.close(fd)
    download(ALL_cCREs_URL, tmpfile_path)

    df = pd.read_table(
        tmpfile_path, 
        header=None,
        names=['Chromosome', 'Start', 'End', 'ID_unknown', 'ID', 'Class'],
    )
    gr = pr.PyRanges(df)
    os.remove(tmpfile_path)

    return gr


AVAILABLE_ORGANS = get_available_organs()


#################################


def read_cCRE_bed(path):
    fd, tmpfile_path = tempfile.mkstemp(suffix='.bed.gz', dir=os.getcwd())
    os.close(fd)
    download(path, tmpfile_path)

    df = pd.read_table(
        tmpfile_path, 
        header=None,
        names=['Chromosome', 'Start', 'End', 'ID', 'Score', 'Strand', 'Start2', 'End2', 'Col', 'Class', 'Completeness'],
    )
    gr = pr.PyRanges(df)
    os.remove(tmpfile_path)

    return gr


def get_cCRE_info_organ(organ):
    assert organ in AVAILABLE_ORGANS, f'"organ" argument must be one of: {AVAILABLE_ORGANS}'

    result = dict()

    http_params = [
        ('type', 'Annotation'),
        ('encyclopedia_version', 'current'),
        ('status', 'released'),
        ('assembly', 'GRCh38'),
        ('annotation_type', 'candidate Cis-Regulatory Elements'),
        ('organism.scientific_name', 'Homo sapiens'),
        ('biosample_ontology.organ_slims', organ),
        ('biosample_ontology.classification', 'primary cell'),
        ('biosample_ontology.classification', 'cell line'),
        ('biosample_ontology.classification', 'in vitro differentiated cells'),
        ('biosample_ontology.classification', 'organoid'),
    ]
    params_string = '&'.join(
        f'{key}={urllib.parse.quote_plus(val)}' 
        for (key, val) in http_params
    )
    query_result = query_ENCODE('search?' + params_string)

    print(f'Found {len(query_result["@graph"])} annotation entries.')
    for item in query_result['@graph']:
        accession = item['accession']
        print(f'Searching for information about {accession}')
        annotation_detail = query_ENCODE(accession)

        biosample_ontology = annotation_detail['biosample_ontology']
        cell_types = biosample_ontology['cell_slims']
        sample_type = biosample_ontology['classification']
        sample_name = biosample_ontology['term_name']

        # get bed file url
        href_candidates = set()
        for file_info in annotation_detail['files']:
            if file_info['annotation_type'] != 'candidate Cis-Regulatory Elements':
                raise Exception(f'"annotation_type" of file entry is not cCRE: {file_info}')
            href = file_info['href']
            if href.endswith('.bed.gz'):
                href_candidates.add(href)
        if len(href_candidates) != 1:
            raise Exception(f'The number of bed file is not 1: {accession}')
        file_url_part = href_candidates.pop()

        sampleinfo = {
            #'accession': accession,
            'cell_types': cell_types,
            'sample_type': sample_type,
            'sample_name': sample_name,
            'file_url': ENCODE_PREFIX + '/' + file_url_part,
        }

        print(f'Downloading cCRE annotation file and loading it as PyRanges object')
        #sampleinfo['bigbed_handle'] = pyBigWig.open(sampleinfo['file_url'])
        sampleinfo['cCRE_gr'] = read_cCRE_bed(sampleinfo['file_url'])
        result[accession] = sampleinfo

    return result
            

def annotate_mutations_with_cCREinfo(coord_list, gr_all_cCRE, organ_cCRE_info, dist=5000):
    # make coord_gr
    chroms = [x[0] for x in coord_list]
    start0s = [x[1] - 1 - dist for x in coord_list]
    end0s = [x[1] + dist for x in coord_list]
    ids = [f'{x[0]}_{x[1]}' for x in coord_list]
    coord_gr = pr.from_dict({'Chromosome': chroms, 'Start': start0s, 'End': end0s, 'Mutation_id': ids})

    # do join operations
    joined_gr_allcCRE = coord_gr.join(gr_all_cCRE)
    joined_grs_organ = dict()
    for accession, sampleinfo in organ_cCRE_info.items():
        joined_grs_organ[accession] = coord_gr.join(sampleinfo['cCRE_gr'])

    # extract results
    result = dict()
    for mutation_id in ids:
        result[mutation_id] = dict()

    for idx, row in joined_gr_allcCRE.df.iterrows():
        cCRE_entry = dict()
        cCRE_entry['start'] = row.Start_b
        cCRE_entry['end'] = row.End_b
        cCRE_entry['chrom'] = row.Chromosome
        cCRE_entry['id'] = row.ID

        cCRE_entry['class'] = dict()
        cCRE_entry['class']['all_samples'] = row.Class

        cCRE_entry['completeness'] = dict()
        cCRE_entry['completeness']['all_samples'] = None

        result[row.Mutation_id][cCRE_entry['id']] = cCRE_entry

    for accession, joined_gr in joined_grs_organ.items():
        sampleinfo = organ_cCRE_info[accession]
        for idx, row in joined_gr.df.iterrows():
            cCRE_entry = result[row.Mutation_id][row.ID]
            #key = (accession, sampleinfo['sample_name'])
            key = accession
            cCRE_entry['class'][key] = row.Class
            cCRE_entry['completeness'][key] = row.Completeness

    return result


def serialize_samples_metadata(organ_cCRE_info):
    data = dict()
    for accession, subdict in organ_cCRE_info.items():
        data[accession] = dict()
        for key in ('cell_types', 'sample_type', 'sample_name'):
            data[accession][key] = subdict[key]

    return json.dumps(data)


#########

def argument_parsing():
    parser = argparse.ArgumentParser(
        description=f'Annotates mutation file with organ-specific ENCODE cCRE information. Only GRCh38 supported. Studies using primary cell, cell line, organoid, or in vitro differentiated cells are searched.'
    )
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--organ', required=True)
    parser.add_argument('--distance', required=False, default=5000, type=int, help=f'cCREs within this distance(unit: base-pair) from input mutations will be searched.')

    args = parser.parse_args()
    return args
        
    
def main():
    args = argument_parsing()

    # downloading ENCODE data files
    print(f'Downloading ENCODE cCRE information files. This may take up to ~10min.')
    print(f'Downloading cCRE bed file for all samples')
    gr_all_cCRE = get_all_cCREs_gr()
    organ_cCRE_info = get_cCRE_info_organ(args.organ)

    # extract mutation coordinates from input file
    print('Beginning annotation')

    infile_df = pd.read_table(args.infile, comment='#', header=None)
    infile_df = infile_df.iloc[:, :2]
    infile_df.columns = ['CHROM', 'POS']

    coord_list = [(row.iloc[0], row.iloc[1]) for idx, row in infile_df.iterrows()]
    mutation_annotation_result = annotate_mutations_with_cCREinfo(coord_list, gr_all_cCRE, organ_cCRE_info, dist=args.distance)
    #mutation_annotation_result_as_string = {key: json.dumps(val) for (key, val) in mutation_annotation_result.items()}
    cCRE_INFO_column = [
        json.dumps(mutation_annotation_result[f'{chrom}_{pos}'])
        for chrom, pos in coord_list
    ]

    # serialize metadata
    samples_metadata_strig = serialize_samples_metadata(organ_cCRE_info)

    # decorate dataframe
    outfile_df = infile_df.assign(cCRE_INFO=cCRE_INFO_column, SAMPLES_METADATA=samples_metadata_strig)

    # write output file
    outfile_df.to_csv(args.outfile, sep='\t', index=False)

    print('All successfully finished')


if __name__ == '__main__':
    main()


