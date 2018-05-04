import os, csv, argparse, sys
import pandas as pd

parser = argparse.ArgumentParser(description='Take csv data from PubMLST and grab what you want out of it.')

parser.add_argument('-infile', metavar='--infile', help='Path to pubMLST info .csv file')
parser.add_argument('-outfile', metavar='--outfile', type=str, help="Output file name")
parser.add_argument('-printids', metavar='--print_ids', type=str, help="Output resulting IDs to specified filename for further analysis. Use .txt or .csv.")
parser.add_argument('-ids', metavar='--ids', nargs='+', type=str, default=None,\
                    help='Give me a list of IDs to include. Mutually exclusive with -idc.')
parser.add_argument('-idc', metavar='--id_csv', nargs='+',type=str, default=None,\
                    help='Give me a comma-separated list of IDs. All on one row is preferable, but if theyre on multiple rows use the -stupid flag. Use .csv or .txt.')
parser.add_argument('-stupid', metavar='It is as it does.', nargs='?', default=None, type=bool,\
                    help='Looks like you put all your IDs on separate rows. That was stupid. You must write "-stupid True"')
parser.add_argument('-rg', metavar='--regions', nargs='+', default=None, type=str, \
                    help='Choose which regions you want to grab samples from.')
parser.add_argument('-c', metavar='--country', nargs='+', default=None, type=str, help='Countries you want to include')
parser.add_argument('-y', metavar='--years', nargs='+', default=None,\
                    help='Which years you want to grab samples from. Use -yr if you have a range of values.')
parser.add_argument('-yr', metavar='--yrange', nargs=2, type=int, default=None,\
                    help='Range of years you want to grab samples from. Provide in this format: "-yr 19XX 20XX"')
parser.add_argument('-sc', metavar='--serotypes_csv', nargs='+', type=str, default=None,\
                    help='Give me a comma-separated list of serotypes. Put them all on the same row.')
parser.add_argument('-s', metavar='--serotypes', nargs='+', type=str, default=None,\
                    help='Type in your serotypes, separated by spaces.')
parser.add_argument('-dia', metavar='--diagnosis', nargs='+', type=str, default=None, help='Diagnoses to include.')
parser.add_argument('-p', metavar='--penicillin', nargs=2, type=float, default=None,\
                    help='Range of penicillin sensitivity values you want. Ex.: -p 0.0 3.0. Decimal point precision enabled.')
parser.add_argument('-e', metavar='--erythromycin', nargs=2, type=float, default=None,\
                    help='Range of erythromycin sensitivity values you want. Ex.: -e 0.0 3.0. Decimal point precsion enabled.')
parser.add_argument('-t', metavar='--tetracycline', nargs=2, type=float, default=None,\
                    help='Range of tetracycline sensitivity values you want. See above examples.')
parser.add_argument('-ch', metavar='--chloramphenicol', nargs=2, type=float, default=None,\
                    help='Range of chloramphenicol sensitivity values you want. See above examples.')
parser.add_argument('-cf', metavar='--cefotaxime', nargs=2, type=float, default=None,\
                    help='Range of cefotaxime sensitivity values you want. See above examples.')
parser.add_argument('-ST', metavar='Desired ST(s)', nargs='+', type=str, default=None, help='Type out the ST types you want.')
parser.add_argument('-MLST', metavar='Desired MLST type', nargs='+', type=str, default=None, help='Type out the MLST types you want.')


args = parser.parse_args()

STUPID = False

infile = args.infile
outfile = args.outfile

ids = args.ids
idc = args.idc

if args.idc is not None:
    if args.ids is not None:
        print("Only provide input to -id or -idc. Not both.")
        sys.exit(420)
    id_csv = args.idc[0]
    ID_CSV = True
    if args.stupid is not None:
        STUPID = True

regions = args.rg
if type(regions) is not list and regions is not None:
    regions = list(regions)

countries = args.c
if type(countries) is not list and countries is not None:
    countries = list(countries)

years = args.y
if type(years) is not list and years is not None:
    years = list(years)

yrange = args.yr
#Insert assert statement

serotypes_csv = args.sc
serotypes = args.s
if serotypes is not None and serotype_csv is not None:
    print("Only provide input to -s or -sc. Not both.")
    sys.exit(666)
if type(serotypes) is not list and serotypes is not None:
    serotypes = list(serotypes)

diagnosis = args.dia
if type(diagnosis) is not list and diagnosis is not None:
    diagnosis = list(diagnosis)

penicillin = args.p
erythromycin = args.e
tetracycline = args.t
chloramphenicol = args.ch
cefotaxime = args.cf

ST = args.ST
if type(ST) is not list and ST is not None:
    ST = list(ST)

MLST = args.MLST
if type(MLST) is not list and MLST is not None:
    MLST = list(MLST)



def getter(field, list, df):
    return df[df[field].isin(list)]

def split_get(field, twolist, df):
    df2 = df.apply(lambda x: x >= twolist[0])
    return df2.apply(lambda x: x <= twolist[1])


try:
    mlst = pd.read_csv(infile, sep=',', low_memory=False)
except:
    print("Error opening your infile. Try again.")

#Let's do all the list-based methods first.
if ids is not None:
    mlst = getter('id', ids, mlst)
if idc is not None and STUPID == False:
    with open(id_csv, 'rb') as infile:
        reader = csv.reader(infile)
        ids = list(reader)
    mlst = getter('id', ids, mlst)
if STUPID:
    id_df = pd.read_csv(id_csv, header=None)
    id_temp = id_csv
    mlst = getter('id', id_df[0], mlst)
if regions is not None:
    mlst = getter('region', regions, mlst)
if countries is not None:
    mlst = getter('country', countries, mlst)
if years is not None:
    mlst = getter('year', years, mlst)
if serotypes is not None:
    mlst = getter('serotype', serotypes, mlst)
if serotypes_csv is not None:
    with open(serotypes_csv, 'wb') as infile:
        reader = csv.reader(infile)
        serotypes = list(reader)
    mlst = getter('serotype', serotypes, mlst)
if ST is not None:
    mlst = getter('ST (MLST)', ST, mlst)
if MLST is not None:
    mlst = getter('clonal_complex (MLST)', MLST, mlst) #MLST is input list; mlst is the dataframe

#Now for the harder ones.
if yrange is not None:
    mlst = split_get('year', yrange, mlst)
if penicillin is not None:
    mlst = split_get('penicillin', penicillin, mlst)
if erythromycin is not None:
    mlst = split_get('erythromycin', erythromycin, mlst)
if tetracycline is not None:
    mlst = split_get('tetracycline', tetracycline, mlst)
if chloramphenicol is not None:
    mlst = split_get('chloramphenicol', chloramphenicol, mlst)
if cefotaxime is not None:
    mlst = split_get('cefotaxime', cefotaxime, mlst)

#Finish up- prepare outfile name

mlst.to_csv(outfile, sep=',')

print('boogie')
