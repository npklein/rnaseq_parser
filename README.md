
INSTALL
-------
python setup.py install

EXAMPLES
--------
from molgenis_api import molgenis
# make a connection
connection = molgenis.Connect_Molgenis('http://localhost:8080', 'admin', 'admin')
# add a row to the entity public_rnaseq_Individuals
connection.add_entity_row('public_rnaseq_Individuals',{'id':'John Doe','age':'26', 'gender':'Male'})
# get the rows from public_rnaseq_Individuals where gender = Male
print connection.get_entity_rows('public_rnaseq_Individuals',[{'field':'gender', 'operator':'EQUALS', 'value':'Male'}])['items']
# update row in public_rnaseqIndivduals where id=John Doe -> set gender to Female
connection.update_entity_row('public_rnaseq_Individuals',[{'field':'id', 'operator':'EQUALS', 'value':'John Doe'}], {'gender':'Female'})









Python command line tool to parse all data from our PublicRNAseq pipeline into a molgenis database

To install parse_rnaseq_output do

sudo python setup.py install

To run do

python run_parser.py --all
