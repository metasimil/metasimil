#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This example runs a similarity analisis on the marcxml collection 2k_a.xml
# 

import metasimil as ms
import codecs
import datetime

# Set up the collection filename
collection = '2k_a.xml'

# Set up the fields that will be analized and thier relative weight in the analysis (weighted geometric mean)
# - The field 'name' is arbitrary (here e.g.: Title, Identifier, Abstract, Date, Creator). Dublincore labels were used.
# - The field 'type' are predefined; appropriated field types will yield better results:
#    - 'wc': adapted for textual fields such as titles or abstracts
#            'wc' means word count: they will be indexed in a dictionary containing the word frequencies,
#             as well as a global word document frequency for the whole collection, required by similarity
#    - 'year' : adapted for date or year fields
#    - 'id' : adapted for unique identifiers fields, such as DOIs, ISBNs, or recids
#    - 'authors' : adapted for people names

fields = { 'Title'      : {'type':'wc'     , 'weight':10  },\
           'Identifier' : {'type':'id'     , 'weight':100 },\
           'Abstract'   : {'type':'wc'     , 'weight':10 },\
           'Date'       : {'type':'year'   , 'weight':2  },\
           'Creator'    : {'type':'authors', 'weight':5  } }

# Loading the data

converterIterator = ms.Converter(collection,format='marcxml')
coll = ms.Collection('ArchiveOuverte', fields=fields)    
coll.converter_load(converterIterator)

# Set ups analysis parameters

# Global similarity output threshold
# Similarity outputs range form 0.0 to 1.0. As a the analyisis of a
# collection consist of huge number of record pairs comaprisons
# introducing a threshold avoids filling up the disk with results.
coll.threshold = 0.75

# The Luhn option will reduce the number of record paris comparisons
# by selecting only tje pairs that ara more likely to contain similar
# records based on a 'wc' index (see above).
# Tiles are generally a good choice for this.

#coll.luhn = False
coll.luhn = 'Title'

# The Luhn speed up will by default do the safe thing : if good candidate
# for similar records cannot be computed for a given record, it is compared
# to all other records of the collection. Exept if the quick and dirty
# option is set. This unsafe option is generally not recommaned.
# It can be usefull for huge collections.
coll.luhn_quick_and_dirty = False
#coll.luhn_quick_and_dirty = True

# Numer of processes, or in practices processors that will be used
# False means just one process.
#coll.processors = False
coll.processors = 3

# Reporting (in a log)
coll.log = 'run.log'
#coll.log = False
coll.detailed_metadata = True
#coll.detailed_metadata = False

# Saving collection to a pickle so it can be reused for a later analyis.
#coll.pickle_dump('archive-ouverte.pickle')

# Analysing and writing results to disk
open('output.csv','w').write('\t'.join(['uuid1','uuid2','similarity','details','recid1','recid2','date1','date2','title1','title2','creator1','creator2'])+'\n')
f = open('output.csv','a')
print('Analyzing collection [%s] ...' % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),) )
for r in coll * coll: # Performs acutal analysis on collection; results are yielded one by one to avoid heavy memory load.
    f.write('\t'.join( [ r['rec1'], r['rec2'],str(r['similarity']), str(r['metadata']),\
                         str(coll.records[r['rec1']].fields['Identifier'].raw), str(coll.records[r['rec2']].fields['Identifier'].raw),\
                         str(coll.records[r['rec1']].fields['Date'].raw), str(coll.records[r['rec2']].fields['Date'].raw),\
                         str(coll.records[r['rec1']].fields['Title'].raw), str(coll.records[r['rec2']].fields['Title'].raw),\
                         str(coll.records[r['rec1']].fields['Creator'].raw), str(coll.records[r['rec2']].fields['Creator'].raw) ] )+'\n')
f.close()
print('Done.')
