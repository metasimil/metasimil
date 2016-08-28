#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This file is part of MetaSimil.
#
# MetaSimil is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MetaSimil is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MarcXimiL.  If not, see <http://www.gnu.org/licenses/>.
#
#
# A copy of the GNU General Public License is located in : ./Licesse/gpl.txt

"""
MetaSimil is a similarity analysis software and Python library.

Date: 2016-08-27
Version: 0.1
Warning: this is an "alpha" release. For production, http://marcximil.sourceforge.net is recomended.

Author: Jan Krause. 
Licence: GNU General Public License v3.0 or later at your option.
Copyrights: Many fucntions have been adapted form MarcXimiL : http://marcximil.sourceforge.net/ (by Alain Borel and Jan Krause, under GNU GPL3.0) 

Strucutre of this library:
The library and the data is organized in collections (class Collection), composed of records (class Records), in turned composed of fields (class Field). An additional Converter class is dedicated to the conversion of metadata formats (CSV, DublinCore, MarcXML) to Python dictionanries. 

Dependencies: Python >= 2.7 or Python >= 3.0.

"""

from math import log, sqrt, log10, exp, e
import copy
import unicodedata
import datetime
from sys import version
from uuid import uuid4
from xml.dom import minidom
from operator import itemgetter
import pickle

if version[0]=='2':
    import codecs

class Field:
    ###############################
    ### text processing helpers ###
    ###############################

    def __remove_accents(self,st):
        """
        This function removes all diacritics.
        """
        if version[0] == '2':
            return unicodedata.normalize("NFKD", unicode(st) ).encode('ascii','ignore')
        else:
            return self.__remove_accents_python3(st)

    def __remove_accents_python3(self,st):
        """
        This function removes lower case diacritics. USED BY PYTHON3 !!!
        """
        accents = { 'áàâä' : 'a' , \
                    'éèêë' : 'e' , \
                    'íìîï' : 'i' , \
                    'óòôö' : 'o' , \
                    'úùûü' : 'u' , \
                    'ýỳŷÿ' : 'y' , \
                    'ç'    : 'c' , \
                    'ñ'    : 'n' }
        output = ''
        for ch in st:
            nch = ch
            for ch_group in accents.keys():
                if ch in ch_group:
                    nch = accents[ch_group]
            output += nch
        return output

    def __remove_tail(self,st):
        """
        This function removes all useless characters at the end of its string argument.
        """
        undesirables = " \t\r\n"
        if len(st) > 0:
            while st[-1] in undesirables:
                st = st[:-1]
                if len(st)==0:
                    break
            return st
        else:
            return ''

    def __remove_head(self,st):
        """
        This function removes all useless characters at the beginning of its string argument.
        """
        undesirables = " \t\r\n"
        if len(st) > 0:
            while st[0] in undesirables:
                st = st[1:]
                if len(st)==0:
                    break
            return st
        else:
            return ''

    def __remove_punctuation(self,st):
        """
        This function removes all useless punctuation signs.
        """
        undesirables = "\'\t\r\n\".;,:!?()"
        for ch in undesirables:
            st = st.replace(ch, '')
        while '  ' in st:
            st = st.replace('  ', ' ')
        return st


    def __normalize_space(self,st):
        """
        This function converts all type of white spaces (tabs, linfeeds, ...) into single spaces.
        """
        undesirables = "\t\r\n"
        for ch in undesirables:
            st = st.replace(ch, ' ')
        return st
        

    def __disambiguate_author(self,st):
        """
        This function tries to disambiguate author names by:
        - normalizing case
        - normalizing diacritics
        - normalizing spaces and punctuation
        - then
            - if there is a coma:
                - it is considerd that the string before the comma is the last name, the rest the firstnames
                - the last name is unaltered, only the inital of the first first name is kept
            - otherwise:
                - takes the last word as the last name, the initial of the first first name
            
        """
        st2 = st.lower()
        if len(st2)>0:
            if ',' in st2:
                ns = st2.split(',')
                fn = self.__normalize_space( self.__remove_head( self.__remove_tail( ns[1] )))
                if len(fn)>0:
                    ifn = fn[0]
                else:
                    ifn = ''
                ln = self.__normalize_space( self.__remove_head( self.__remove_tail( ns[0] )))
            else:
                stx = self.__normalize_space( self.__remove_head( self.__remove_tail( st2 )))
                nx = stx.split(' ')
                if len(nx[0])>0:
                    ifn = nx[0][0]
                else:
                    ifn = ''
                if len(nx[-1])>0:
                    ln = nx[-1]
                else:
                    ln = ''
            return self.__remove_accents(ln + ' ' + ifn)
        else:
            return None

    def __text2flatlist(self, st):
        """
        This function transforms a text into a list:
        - punctuation is removed, including parentheses and apostrophes (square brackets are kept)
        - double spaces, tabs, linefeeds are removed
        - everything is put to lower case
        - accents are removed
        """
        if st == None or st == '' or st == ' ':
            return None
        stl = st.lower()
        stm = str( self.__remove_accents( self.__remove_punctuation( self.__normalize_space( self.__remove_tail( self.__remove_head(stl ))))) )
        stn = stm.split(' ')
        if '' in stn:
            stn.remove('')
        return stn

    def __frequencies(self, lst, options = None):
        """
        This function returns a dictionary of the frequencies of items contained in a list.
        For exemple the list: ['a', 'a', 'a', 'b', 'c'] would yield {'a':3, 'b':1, 'c':1 }
        """
        if lst == None:
            return None
        items = set(lst)
        output = {}
        for element in items:
            output[element] = lst.count(element)
        return output

    ############################
    ### similariy algorithms ###
    ############################


    def __get_docfrequency(self, parentcoll, fieldname, t1):
        if t1 in parentcoll.global_wc[fieldname].keys():
            return parentcoll.global_wc[fieldname][t1]['df']
        else:
            return None
        
    def __update_docfrequency(self, parentcoll, fieldname, t1, count):
        parentcoll.global_wc[fieldname][t1]['df'] += count

    def __get_nrecords(self, parentcoll):
        return parentcoll.size

    #####################################
    ### records comparison functions  ###
    #####################################

    def __ntfnidf_vectorcosine_comp_obj__wc(self,obj1, obj2):
        """
        Calculates the vector similarity of two word count dictionaries
        using the cosine definition
        and Salton & Turtle's ntf*nidf weight coefficients
        """
        
        if (obj1.wc == None) or (obj2.wc == None) or (len(obj1.wc) == 0) or (len(obj2.wc) == 0) :
            return None

        maxtf1=max(obj1.wc.values())*1.0
        w1={}
        maxtf2=max(obj2.wc.values())*1.0
        w2={}
        n=self.__get_nrecords(self.parentcoll)

        sqnorm1=0.0
        for t1 in obj1.wc:
            # if comparing fields of same coll use coll docfrequency
            # if comapring fields in different colls, use sum of docfrequency
            if obj1.parentcoll == obj2.parentcoll:             
                df = obj1.__get_docfrequency(self.parentcoll,self.fieldname,t1)
            else:
                df1 = obj1.__get_docfrequency(obj1.parentcoll,self.fieldname,t1)
                df2 = obj2.__get_docfrequency(obj2.parentcoll,self.fieldname,t1)
                if (df1 is not None) and (df2 is not None):
                    df = max([df1, df2])
                elif df1 is not None:
                    df = df1
                elif df2 is not None:
                    df = df2
                else:
                    df = None
            try: #FIXME-float-div-0
                w1[t1] = obj1.wc[t1]/maxtf1 * log(n/df)/log(n)
            except:
                w1[t1] = 0
            sqnorm1 += w1[t1]*w1[t1]

        sqnorm2=0.0
        dot_product=0.0
        for t2 in obj2.wc:
            try: #FIXME-float-div-0
                w2[t2] = obj2.wc[t2]/maxtf2 * log(n/obj1.__get_docfrequency(self.parentcoll,self.fieldname,t2))/log(n)
            except:
                w2[t2] = 0
            sqnorm2 += w2[t2]*w2[t2]
            try:
                dot_product += w2[t2]*w1[t2] 
            except KeyError:
                pass
        
        if (sqnorm1 * sqnorm2) == 0:
            return 0
        else:
            return dot_product/sqrt(sqnorm1*sqnorm2)

    ###################################
    ### Fields comparison functions ###
    ###################################
    
    def __years_comp_markximil__raw(self, obj1, obj2, options = None):
        """
        This function tries to extract and then compare the years in format YYYY
        between two fields.
        - if the years are the same 1.0 is returned
        - for each year of difference 0.1 is substracted to the 1.0 similarity
        - if anything fails None is returned
        """
        if not( isinstance(obj1.raw, list) and isinstance(obj2.raw, list) ):
            return None #FIXME:None
        
        f1s = str(obj1.raw[0])
        f2s = str(obj2.raw[0])
        try:
            nd1 = int(f1s[:4])
            nd2 = int(f2s[:4])
        except:
            return None #FIXME:None
            
        if nd1 == nd2:
            return 1.0
        else:
            return 0.5 ** ( abs(nd1-nd2) )

    def __years_comp__raw(self, obj1, obj2, options = None):
        """
        This function tries to extract and then compare the years in format YYYY
        between two fields:
        - if the years are the same 1.0 is returned
        - for each year of difference 0.1 is substracted to the 1.0 similarity
        - if anything fails None is returned
        """
        if not( isinstance(obj1.raw, list) and isinstance(obj2.raw, list) ):
            return None #FIXME:None
        
        f1s = str(obj1.raw[0])
        f2s = str(obj2.raw[0])
        try:
            nd1 = int(f1s[:4])
            nd2 = int(f2s[:4])
        except:
            return None #FIXME:None
            
        if nd1 == nd2:
            return 1.0
        elif abs(nd1-nd2) < 10:
            return 1.0 - ( 0.1 * abs(nd1-nd2) )
        else:
            return 0.0

        
    def __identifiers_comp_marximil__raw(self, obj1, obj2, options = None):
        """
        This function compares two fields containing identifiers, like:
        - standard digital indentifiers: DOIs, PMIDs, ISSNs, ISBNs, etc.
        - record numbers
        - OAIids

        It is case insensitive.

        The following substrings are eliminated: 'doi:' 'pmid:' 'issn:' 'isbn:' '-'
        
        The output is:
        - 1.0 if the two identifiers are identical
        - 0.0 otherwise
        """
        unwanted = ['doi:', 'pmid:', 'issn:', 'isbn:', '-', 'oai:']

        f1 = obj1.raw
        f2 = obj2.raw
        
        if (f1 == None) or (f2 == None):
            return None
        else:
            id1 = self.__remove_head( self.__remove_tail( f1 )).lower()
            id2 = self.__remove_head( self.__remove_tail( f2 )).lower()
            for x in unwanted:
                id1 = id1.replace(x, '')
                id2 = id2.replace(x, '')
            if id1 == id2 :
                return 1.0
            else:
                return 0.0


    def __multiidentifiers_comp__raw(self, obj1, obj2, options = None):
        """
        This function compares two sets of identifiers and looks for any shared value.
        Since some identifier conventions are case-insensitive, we convert everything to lowercase.
        We exclude a prefix of 5 chars or less, followed by a colon: it could be
        just an indication of the identifier type (such as DOI:, ArXiv:, etc), it's easier to neglect it.
        """
        
        f1 = obj1.raw
        f2 = obj2.raw

        if (f1 == None) or (f2 == None):
            return None

        set1 = set()
        for x in f1:
            if x is not None:
                colon = x.find(':')
                if colon < 0 or colon > 6:
                    set1 |= set([self.__remove_head( self.__remove_tail(x.lower()))])
                else:
                    set1 |= set([self.__remove_head( self.__remove_tail(x.lower().split(':')[1]))])

        set2 = set()
        for x in f2:
            if x is not None:
                colon = x.find(':')
                if colon < 0 or colon > 5:
                    set2 |= set([self.__remove_head( self.__remove_tail(x.lower()))])
                else:
                    set2 |= set([self.__remove_head( self.__remove_tail(x.lower().split(':')[1]))])

        if len(set1 & set2) > 0:
            return 1.0
        else:
            return None

    def __authors_comp__raw(self, obj1, obj2, options = None):
        """
        Compare authors lists. Returns percentage of similar authors,
        after disambiguation.
        """
        f1 =obj1.raw
        f2 =obj2.raw
        
        if (f1 == None) or (f2 == None):
            return None
        else:
            fb1 = map(self.__disambiguate_author, f1)
            fb2 = map(self.__disambiguate_author, f2)
            it1 = set(fb1)
            it2 = set(fb2)
            percentage = float( len( it1 & it2 ) ) / len( it1 | it2 ) 
            return (exp(percentage**(1.8)) -1) / (e-1) #1.7183
        
    ###################
    ### constructor ###
    ###################
    
    def __init__(self, parentcoll, fieldname, value):
        # World count fields need to be text, list are concatenated.
        if parentcoll.fields[fieldname]['type'] == 'wc':
            if value is not None:
                self.raw = ' '.join(value)
            else:
                self.raw = None
        else:
            self.raw = value
        self.fieldname = fieldname
        self.parentcoll = parentcoll
        if parentcoll.fields[fieldname]['type'] == 'wc':
            self.wc = self.__frequencies( self.__text2flatlist(self.raw) )
        else:
            self.wc = None

    #######################
    ### overiding stuff ###
    #######################

    def __mul__(self, other):
        #print(self.parentcoll.fields[self.fieldname]['type'])
        if self.parentcoll.fields[self.fieldname]['type'] == 'wc':
            return self.__ntfnidf_vectorcosine_comp_obj__wc(self, other)
        elif self.parentcoll.fields[self.fieldname]['type'] == 'year':
            return self.__years_comp__raw(self, other)
        elif self.parentcoll.fields[self.fieldname]['type'] == 'id':
            return self.__multiidentifiers_comp__raw(self, other)
        elif self.parentcoll.fields[self.fieldname]['type'] == 'authors':
            return self.__authors_comp__raw(self, other)
        else:
            return None # FIXME should be None, shouldn't it?

    def __repr__(self):
        return self.raw

    def __str__(self):
        return self.raw

class Record:

    ###################
    ### constructor ###
    ###################
    
    def __init__(self, recid, parentcoll, record=None):
        self.parentcoll = parentcoll
        self.recid = recid
        if record is None:
            self.fields = copy.deepcopy(parentcoll.record_template)
        else:
            self.fields = record
        parentcoll.size = len(parentcoll.records)
   
    def set_field(self, fieldname, value):
        """ If field does not exist, the field is created. Otherwise it is updated."""
        # if the field does not exist it is set to None -> creating
        if self.fields[fieldname] == None:
            self.fields[fieldname] = Field(self.parentcoll, fieldname, value)
            # if fieldtype is word count (wc), global wc index is updated
            if self.parentcoll.fields[fieldname]['type'] == 'wc':
                if self.fields[fieldname].wc is not None: #FIXME
                    for k in self.fields[fieldname].wc.keys():
                        ks = self.parentcoll.global_wc[fieldname].keys()
                        if not k in ks:
                            self.parentcoll.global_wc[fieldname][k] = {'df':self.fields[fieldname].wc[k], 'docs':[self.recid]}
                        else:
                            self.parentcoll.global_wc[fieldname][k]['df'] += self.fields[fieldname].wc[k]
                            self.parentcoll.global_wc[fieldname][k]['docs'].append(self.recid)
        # else updating
        else:
            if self.parentcoll.fields[fieldname]['type'] == 'wc':
                old_wc = self.fields[fieldname].wc
            self.fields[fieldname] = Field(self.parentcoll, fieldname, value)
            # if fieldtype is word count (wc), global wc index is updated
            if self.parentcoll.fields[fieldname]['type'] == 'wc':
                for k in self.fields[fieldname].wc.keys():
                    ks = self.parentcoll.global_wc[fieldname].keys()
                    if not k in ks:
                        self.parentcoll.global_wc[fieldname][k] = {'df':self.fields[fieldname].wc[k], 'docs':[self.recid]}
                    else:
                        self.parentcoll.global_wc[fieldname][k]['df'] += self.fields[fieldname].wc[k]
                        if k in old_wc.keys():
                            self.parentcoll.global_wc[fieldname][k]['df'] -= old_wc[k]['df']

    def remove_field(self, fieldname):
        """ Deletes a field, updates global-wc if needed."""
        if fieldname in self.fields.keys():
            if self.parentcoll.fields[fieldname]['type'] == 'wc':
                for k in self.fields[fieldname].wc.keys():
                    self.parentcoll.global_wc[fieldname][k] -= self.fields[fieldname].wc[k]
            del(self.fields[fieldname])

                            
    #######################
    ### overiding stuff ###
    #######################

    def __mul__(self, other):
        # performing weighted geometric mean
        eps = 0.00000000001
        sf = 0
        total_weight = 0
        if not self.parentcoll.detailed_metadata:
            dmd = False
            metadata = None
        else:
            dmd = True
            metadata = {}
        for k in self.fields.keys():
            if (self.fields[k] is not None) and (other.fields[k] is not None):
                field_sim = self.fields[k] * other.fields[k]
                if dmd:
                    metadata[k] = field_sim
                #print(field_sim)
                if field_sim is not None:
                    total_weight += self.parentcoll.fields[k]['weight']
                    sf += self.parentcoll.fields[k]['weight'] * log(field_sim + eps)
        if total_weight > 0:
            return {'similarity':exp( sf / total_weight ) - eps , 'metadata':metadata}
        else:
            return {'similarity':None, 'metadata':metadata}

    def __repr__(self):
        return self.recid

    def __str__(self):
        return  self.recid

    
class Collection:

    ###################
    ### constructor ###
    ###################
    
    def __init__(self, name, fields = { 'title':{'type':'wc', 'coeff':1} } ):

        # Name of the collection
        self.name = name        
        self.fields = fields
        self.threshold = 0.005
        self.luhn = False
        self.luhn_quick_and_dirty = False
        self.detailed_metadata = False
        self.processors = False
        self.log = False

        # Initializing records, a dict with recids as keys
        self.records = {} 
        self.size = len(self.records)

        # Initializing global word count index
        self.global_wc = {}
        for k in fields.keys():
            if self.fields[k]['type']=='wc':
                self.global_wc[k] = {}

        # Initializing record template
        self.record_template = {}
        for k in self.fields.keys():
            self.record_template[k] = None

    def set_record(self, recid, record = None):
        """ If record is given (not None), appropriated fields will be loaded. If record exists it is updated"""
        self.records[recid] = Record(recid, self, record=None)
        if record is not None:
            # if record is not None addin fileds common to
            # given record and field description
            #print(record)
            kr = record.keys()
            kd = self.fields.keys()
            common_keys = list( set(kr) & set(kd) )
            for k in common_keys:
                self.records[recid].set_field(k, record[k])

    def remove_record(self, recid):
        """ Deletes a record, after deletion of all its fields (global-wc are updated with that operation)"""
        ks = self.records[recid].fields.keys()
        for fieldname in list(ks):
            self.records[recid].remove_field(fieldname)        
        del(self.records[recid])
                
    def compare_records(self):
        """ Comapres together all records within a collection. """
        coll_keys = list(self.records.keys())
        for i in range(len(coll_keys)):
            for j in range(i+1, len(coll_keys)):
                result = self.records[coll_keys[i]] * self.records[coll_keys[j] ]
                if result['similarity'] >= self.threshold:
                    yield { 'rec1': self.records[coll_keys[i]].recid,\
                            'rec2': self.records[coll_keys[j]].recid,\
                            'similarity': result['similarity'],\
                            'metadata': result['metadata']}

    def compare_records__multipross_yield(self):
        """ Comapres together all records within a collection. Multiprocessing. """

        from multiprocessing import Pool
        #from functools import partial
        #compare_pair_partial = partial(compare_pair,self)

        def thresholder(result):
            if result['similarity'] > self.threshold:
                return True
            else:
                return False

        recid_pairs = []
            
        coll_keys = list(self.records.keys())
        
        for i in range(len(coll_keys)):
            for j in range(i+1, len(coll_keys)):
                recid_pairs.append((coll_keys[i],coll_keys[j]))

        # yielding results by chunks
    
        pool = Pool(processes=self.processors)
        n_pairs = len(recid_pairs)
        chunk_size = 1000000

        if n_pairs <= chunk_size:

            chunk_pairs = []
            for x in recid_pairs:
                # FIXME : possible improvement
                # avoid passing whole records, but jeust keys : memory usage
                # chunk_pairs.append( (x[0], x[1]) )
                chunk_pairs.append( (self.records[x[0]], self.records[x[1]]) )

            for x in pool.map(compare_pair, chunk_pairs):
                if thresholder(x):
                    yield x
        else:

            n_chunks = n_pairs // chunk_size
            
            for chunk in range(n_chunks):
                # extrcating recid pairs for this chund
                if not chunk == (n_chunks-1):
                    chunk_recid_pairs = recid_pairs[chunk*chunk_size:(chunk+1)*chunk_size ]
                else:
                    chunk_recid_pairs = recid_pairs[chunk*chunk_size:]

                # getting records form recids
                chunk_pairs = []
                for x in chunk_recid_pairs:
                    #chunk_pairs.append( (x[0], x[1]) )
                    chunk_pairs.append( (self.records[x[0]], self.records[x[1]]) )

                #for x in pool.map(compare_pair_partial, chunk_pairs):
                for x in pool.map(compare_pair, chunk_pairs):
                    if thresholder(x):
                        yield x

    def compare_records__multipross(self):
        """ Comapres together all records within a collection. Multiprocessing. 
            All results are yielded in one batch, witch may lead to a lot
            of memory usage. compare_records__multipross fixes that problem."""

        def thresholder(result):
            if result['similarity'] > self.threshold:
                return True
            else:
                return False

        from multiprocessing import Pool
        #from functools import partial
        #compare_pair_partial = partial(compare_pair,self)
        
        pairs = []
        
        coll_keys = list(self.records.keys())
        for i in range(len(coll_keys)):
            for j in range(i+1, len(coll_keys)):
                pairs.append((self.records[coll_keys[i]],self.records[coll_keys[j]]))
                #pairs.append((coll_keys[i],coll_keys[j]))

        with Pool(processes=self.processors) as pool:
            #for x in pool.map(compare_pair_partial, pairs):
            for x in pool.map(compare_pair, pairs):
                if thresholder(x):
                    yield x
    
    def compare_records_two_collections(self, other):
        """ Comapres together all records of the first collectio to all records of the second. """
        coll1_keys = list(self.records.keys())
        coll2_keys = list(other.records.keys())
        threshold = max([self.threshold, other.threshold])
        for i in range(len(coll1_keys)):
            for j in range(len(coll2_keys)):
                result = self.records[coll1_keys[i]] * other.records[coll2_keys[j] ]
                if result['similarity'] >= threshold:
                    yield { 'rec1': self.records[coll1_keys[i]].recid,\
                            'rec2': other.records[coll2_keys[j]].recid,\
                            'similarity': result['similarity'],\
                            'metadata': result['metadata']}

    def compare_records__luhn(self):
        """ Comapres together all records within a collection, 
            only if titles have common significant words. 
            This requires that a field named "Title" of type wc exists"""

        luhn_field = self.luhn
        
        coll_keys = list(self.records.keys())
        for i in range(len(coll_keys)):

            # luhn usefull words
            # frist getting docfrequency for each word in the lunh_field
            # it is sotred in the coll global index: print(self.global_wc[luhn_field])
            title_index = {}

            if self.records[coll_keys[i]].fields[luhn_field].wc is not None:    
                for word in self.records[coll_keys[i]].fields[luhn_field].wc.keys():
                    title_index[word] = self.global_wc[luhn_field][word]
            else:
                 title_index = {}  

            if version[0]=='2':
                useful_words = [x for x in sorted(title_index.items(), key=itemgetter(1)) if x[1]['df'] > 1]
            else:
                useful_words = [x for x in sorted(title_index.items(), key=lambda z: z[1]['df']) if x[1]['df'] > 1]

            if self.luhn_quick_and_dirty:
                if (len(useful_words) == 0):
                    to_check = []
                # if there is 1 word we use it 
                elif (len(useful_words) == 1):
                    # the records to compare to all other records containing it
                    list1 = set(title_index[useful_words[0][0]]['docs']) #list(set(title_index[useful_words[0][0]]['docs']))
                    to_check = [x for x in list1 if coll_keys.index(x) > i]
                    del(list1)
                else:
                    # if there are 2 words or more, we identify the two words closest to median df
                    below_median = (len(useful_words)-1)//2
                    # the records to compare must contain at least one of these words
                    set1=set(title_index[useful_words[below_median][0]]['docs'])
                    set2=set(title_index[useful_words[below_median+1][0]]['docs'])
                    orlist = set1 | set2 #list(set1 | set2)
                    del(set1)
                    del(set2)
                    to_check = [x for x in orlist if coll_keys.index(x) > i]
                    del(orlist)
            else:
                if (len(useful_words) < 2):
                    to_check_range = range( i+1, len(coll_keys) )
                    to_check = [coll_keys[i] for i in to_check_range]
                else:
                    # if there are 2 words or more, we identify the words closest to the median df
                    below_median = (len(useful_words)-1)//2
                    # the records to compare must contain at least one of these words
                    set1=set(title_index[useful_words[below_median][0]]['docs'])
                    set2=set(title_index[useful_words[below_median+1][0]]['docs'])
                    orlist = set1 | set2 #list(set1 | set2)
                    del(set1)
                    del(set2)
                    to_check = [x for x in orlist if coll_keys.index(x) > i]
                    del(orlist)
                
            # # if there is only 0 or 1 title word, we compare with the whole collection
            # if (len(useful_words) < 2):
            #     to_check_range = range( i+1, len(coll_keys) )
            #     to_check = [coll_keys[i] for i in to_check_range]
            # else:
            # #elif False:
            #     # if there are 2 words or more, we identify the words with a median df
            #     below_median = (len(useful_words)-1)//2
            #     #print (len(useful_words),below_median)
            #     #middle_words = [useful_words[below_median][0],useful_words[below_median+1][0]]
            #     #print(ri, middle_words)
            #     # the records to compare must contain at least one of these words
            #     set1=set(title_index[useful_words[below_median][0]]['docs'])
            #     set2=set(title_index[useful_words[below_median+1][0]]['docs'])
            #     #print(set1,set2)
            #     orlist = list(set1 | set2)
            #     #print('\n'.join(orlist))
            #     to_check = [x for x in orlist if coll_keys.index(x) > i]
            #     # FIXME
            #     #print(to_check)
            #     #print(coll_keys)

            for kk in to_check:
                #print(coll_keys[i], kk)
                result = self.records[coll_keys[i]] * self.records[kk]
                if result['similarity'] is not None:
                    if result['similarity'] >= self.threshold:
                        yield { 'rec1': self.records[coll_keys[i]].recid,\
                                'rec2': kk,\
                                'similarity': result['similarity'],\
                                'metadata': result['metadata']}

    def compare_records__luhn__multipross_yield(self):
        """ Comapres together all records within a collection, 
            only if titles have common significant words. 
            This requires that a field named "Title" of type wc exists.
            The records to compare are selected according 2 criteria:
            - the record title must contain one or two words of the reference record title,
             chosen for their median or near-median document frequency (chosen for discriminating power, as in Luhn's Law)
             - the record number must be higher than the reference record's
            This function is compatible with caching, but can also be used without it.

        """

        luhn_field = self.luhn
        
        from multiprocessing import Pool
        from functools import partial
        #compare_pair_partial = partial(compare_pair,self)

        def thresholder(result):
            if result['similarity'] > self.threshold:
                return True
            else:
                return False

        recid_pairs = []
        
        coll_keys = list(self.records.keys())
        
        if self.log:
            write_message(self.log, 'Luhn speedup starting.')
            write_message(self.log, 'Number of records: '+str(len(coll_keys)))
            write_message(self.log, 'Computing pairs using Luhn useful words...')

        pool = Pool(processes=self.processors)

        chunk_size = 2000#0

        luhn_evals = range(len(coll_keys))
        luhn_evals_chunks = []
        for i in range(0, len(luhn_evals), chunk_size):
            luhn_evals_chunks.append( luhn_evals[i:i + chunk_size] )
        del(luhn_evals)
            
        for chunk in luhn_evals_chunks:
            
            if self.log:
                write_message(self.log, '    Computing luhn pairs through colection %s / %s' % (i, len(coll_keys),) )

            func = partial(luhn_eval,self)
            pair_lists = pool.map(func, chunk)
            for chk in pair_lists:
                for pair in chk:        
                    recid_pairs.append(pair)

        chunk_size = 200000#0
        n_pairs = len(recid_pairs)

        if n_pairs <= chunk_size:

            if self.log:
                write_message(self.log,'    Processing all pairs in one go...')
            
            chunk_pairs = []
            for x in recid_pairs:
                chunk_pairs.append( (self.records[x[0]], self.records[x[1]]) )
                #chunk_pairs.append( (x[0], x[1]) )

            #for x in pool.map(compare_pair_partial, chunk_pairs):
            for x in pool.map(compare_pair, chunk_pairs):
                if thresholder(x):
                    yield x
        else:

            if self.log:
                write_message(self.log, '    Processing all pairs in chunks...')
            
            n_chunks = n_pairs // chunk_size
            pool = Pool(processes=self.processors)

            if self.log:
                write_message(self.log, '   Number of chunks: %s' % (n_chunks,))

            for chunk in range(n_chunks):
                if self.log:
                    write_message(self.log, '        Processing similarity analysis on chunk %s / %s' % (chunk+1, n_chunks))
                # extrcating recid pairs for this chunk
                if not chunk == (n_chunks-1):
                    chunk_recid_pairs = recid_pairs[chunk*chunk_size:(chunk+1)*chunk_size ]
                else:
                    chunk_recid_pairs = recid_pairs[chunk*chunk_size:]

                # getting records form recids
                chunk_pairs = []
                for x in chunk_recid_pairs:
                    #chunk_pairs.append( (x[0], x[1]) )
                    chunk_pairs.append( (self.records[x[0]], self.records[x[1]]) )

                #for x in pool.map(compare_pair_partial, chunk_pairs):
                for x in pool.map(compare_pair, chunk_pairs):
                    if thresholder(x):
                        yield x

                        
    #############################
    ### loading and exporting ###
    #############################
    
    def converter_load(self, converter):
        """ Converter is an object of the Converter class """
        if self.log:
            write_message(self.log, 'Converting collection records (converter_load).')
        for record in converter.yield_records():
            recid = str(uuid4())
            self.set_record(recid, record)

    def pickle_dump(self, filename):
        """ Dumps Collection In File Pickle"""
        if self.log:
            write_message(self.log, 'Dumping collection (pickle_dump). Filename: ' + filename)
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def pickle_load(self, filename):
        """ Loads Collection From Pickle File"""
        if self.log:
            write_message(self.log, 'Loading collection (pickle_load). Filename: ' + filename)
        with open(filename, 'rb') as f:
            self = pickle.load(f)


    #######################
    ### overiding stuff ###
    #######################

    def __mul__(self, other):
        if self.log:
            write_message(self.log, 'Comparing collection(s) (Collection.__mul__): %s AND %s' % (self.name, other.name,))
            write_message(self.log, '    self.processors =  %s' % (self.processors,))
            write_message(self.log, '    self.luhn =  %s' % (self.processors,))   
        if self == other:
            if self.processors and self.luhn:
                return self.compare_records__luhn__multipross_yield()
            elif self.processors:
                return self.compare_records__multipross_yield()
            elif self.luhn:
                return self.compare_records__luhn()
            else:
                return self.compare_records()
        else:
            return self.compare_records_two_collections(other)
        
    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

class Converter:
    """ Conversion of metadata formats in srtigs (MarcXML, JSON, DublinCore) in python dictionaries """
    def __init__(self, filename, format = 'dublincore', options = None):
        self.filename = filename
        self.format =  format
        self.options = options
        self.result = {}

    ####################
    ### File parsing ###
    ####################

    def yield_records(self):
        if self.format == 'marcxml':
            if version[0]=='2':
                with codecs.open( self.filename, "r", "utf-8" ) as f:
                    record = ''
                    aggregating = False
                    for line in f:
                        if '<record' in line:
                            record += line
                            aggregating = True
                        elif '</record' in line:
                            aggregating = False
                            record += line
                            recout = record
                            record = ''
                            yield self.convert(recout)
                        elif aggregating:
                            record += line
            else:            
                with open(self.filename) as f:
                    record = ''
                    aggregating = False
                    for line in f:
                        if '<record' in line:
                            record += line
                            aggregating = True
                        elif '</record' in line:
                            aggregating = False
                            record += line
                            recout = record
                            record = ''
                            yield self.convert(recout)
                        elif aggregating:
                            record += line

        elif self.format == 'dublincore':
            if version[0]=='2':
                with codecs.open( self.filename, "r", "utf-8" ) as f:
                    record = ''
                    aggregating = False
                    for line in f:
                        line2 = line.replace('<dc:','<').replace('</dc:','</')
                        if '<resource' in line2:
                            record += line2
                            aggregating = True
                        elif '</resource' in line2:
                            aggregating = False
                            record += line2
                            recout = record
                            record = ''
                            yield self.convert(recout)
                        elif aggregating:
                            record += line2
            else:
                with open(self.filename) as f:
                    record = ''
                    aggregating = False
                    for line in f:
                        line2 = line.replace('<dc:','<').replace('</dc:','</')
                        if '<resource' in line2:
                            record += line2
                            aggregating = True
                        elif '</resource' in line2:
                            aggregating = False
                            record += line2
                            recout = record
                            record = ''
                            yield self.convert(recout)
                        elif aggregating:
                            record += line2


        elif self.format == 'csv':
            if version[0]=='2':
                with codecs.open( self.filename, "r", "utf-8" ) as f:
                    header = True                    
                    for line in f:
                        if not header:
                            yield self.convert(line)
                        if header:
                            header = False
            else:
                with open(self.filename) as f:
                    header = True                    
                    for line in f:
                        if not header:
                            yield self.convert(line)
                        if header:
                            header = False

    #######################
    ### MARCXML methods ###
    #######################

    def marc_parse_controlfield(self, marc, record):
        """
        This function extracts the content of a control field
        - record is an xmlstring contianing a record
        - marc a string representing a MARC control field: only 3 digits!
        This function returns a unicode string : the field content
        """
        parsed_field = None

        xrec = minidom.parseString(record.encode("utf-8"))
        controlfields = xrec.getElementsByTagName('controlfield')

        for controlfield in controlfields:
            if controlfield.attributes['tag'].value == marc: 
                parsed_field = controlfield.childNodes[0].data

        return [parsed_field]


    def marc_parse_multi(self, marclist, record):
        """
        This function extracts multiple subfields together (ex: all 100__a and all 700__a ).
        - record is an xmlstring contianing a record
        - marclist a list representing all MARC fields that we want to extract
        This function returns a list of unicode strings : the fields contents.
        These fields are returned in the same order as the fields are specified.
        """
        output = []
        #print(record+'\n\n\n')

        xrec = minidom.parseString(record.encode("utf-8"))
        datafields = xrec.getElementsByTagName('datafield')

        for datafield in datafields:
            for marc in marclist:
                try:
                    if datafield.attributes['tag'].value == marc[:3] and \
                            datafield.attributes['ind1'].value == marc[3] and \
                            datafield.attributes['ind2'].value == marc[4]:
                        subfields = datafield.getElementsByTagName('subfield')
                        for subfield in subfields:
                            if subfield.attributes['code'].value == marc[-1]:
                                output.append(subfield.childNodes[0].data)
                except:
                    pass
        if output==[]:
            return None
        else:
            return output
            
    def dublincore_parser(self, record):
        """
        """
        output = {}

        xrec = minidom.parseString(record.encode("utf-8"))

        fieldnames = ['Creator', 'Identifier', 'Date', 'Title', 'Description', 'Subject', 'Contributor,'\
                      'Type', 'Format', 'Source', 'Language', 'Relation', 'Coverage', 'Rights']

        for fldname in fieldnames:
            fields = xrec.getElementsByTagName(fldname.lower())
            output[fldname] = []
            for field in fields:
                output[fldname].append(field.childNodes[0].data)
            if output[fldname] == []:
                del(output[fldname])
        if output=={}:
            return None
        else:
            return output

    def csv_parse_single(self, column, record, sep='\t', delim='"'):
        """
        Parses a filed (column) from a csv record
        """
        parsed_field = None

        rs = record.split(sep)

        if len(rs) >= column:
            parsed_field = rs[column]
            if len(parsed_field) == 0:
                parsed_field = None
            else:
                if parsed_field[0] == delim:
                    parsed_field = parsed_field[1:]
                    if len(parsed_field) > 0:
                        if parsed_field[-1] == delim:
                            parsed_field = parsed_field[:-1]

        return [parsed_field]


    def convert(self, text_metadata):
        if self.format == 'marcxml':
            if self.options == None:
                options = [ {'Source':['245  a'], 'Python':'Title', 'Method':self.marc_parse_multi},
                            {'Source':['260  c'], 'Python':'Date', 'Method':self.marc_parse_multi},
                            {'Source':['100 a', '700  a'], 'Python':'Creator', 'Method':self.marc_parse_multi},
                            {'Source':'001', 'Python':'Identifier', 'Method':self.marc_parse_controlfield} ]
            else:
                options = self.options
            for o in options:
                self.result[o['Python']] = o['Method'](o['Source'], text_metadata)
        elif self.format == 'dublincore':
            self.result = self.dublincore_parser(text_metadata)
        if self.format == 'csv':
            if self.options == None:
                options = [ {'Source':1, 'Python':'Title', 'Method':self.csv_parse_single},
                            {'Source':2, 'Python':'Date', 'Method':self.csv_parse_single},
                            {'Source':3, 'Python':'Creator', 'Method':self.csv_parse_single},
                            {'Source':0, 'Python':'Identifier', 'Method':self.csv_parse_single} ]
            else:
                options = self.options
            for o in options:
                self.result[o['Python']] = o['Method'](o['Source'], text_metadata)
        else:
            pass
        return self.result
    
#####################################
### Multiprocessing Functions     ###
### musst be defiend at top level ###
#####################################

def compare_pair(pair):
    """ Called by mupltiprocessin comparisom functions """
    result = pair[0] * pair[1]
    return  { 'rec1': pair[0].recid,\
              'rec2': pair[1].recid,\
              'similarity': result['similarity'],\
              'metadata': result['metadata']}

# def compare_pair(coll, pair):
#     """ Called by mupltiprocessin comparisom functions """
#     result = coll.records[pair[0]] * coll.records[pair[1]]
#     if True:
#     #if result['similarity'] >= self.threshold:
#         #print(result['similarity'])
#         return  { 'rec1': pair[0],\
#                   'rec2': pair[1],\
#                   'similarity': result['similarity'],\
#                   'metadata': result['metadata']}
    
def luhn_eval(coll,i):
    """ Compute Luhn Pairs for record i. Goal: Reducing the number of comarisons for speedup.
        This function if only used by the multiprocessing variant of Luhn comparisions.
    """

    # luhn usefull words
   
    # frist getting docfrequency for each word in the lunh_field
    # it is sotred in the coll global index: print(self.global_wc[luhn_field])

    recid_pairs = []
    title_index = {}

    coll_keys = list(coll.records.keys())
    luhn_field = coll.luhn

    if coll.records[coll_keys[i]].fields[luhn_field].wc is not None:    
        for word in coll.records[coll_keys[i]].fields[luhn_field].wc.keys():
            title_index[word] = coll.global_wc[luhn_field][word]
    else:
         title_index = {}  

    if version[0]=='2':
        useful_words = [x for x in sorted(title_index.items(), key=itemgetter(1)) if x[1]['df'] > 1]
    else:
        useful_words = [x for x in sorted(title_index.items(), key=lambda z: z[1]['df']) if x[1]['df'] > 1]

    # MarcXimiL: if there is only 0 or 1 title word, we compare with the whole collection
    # Metasim QuickAndDirty option: if there is 0 word we compare with the whole collection
    if coll.luhn_quick_and_dirty:
        if (len(useful_words) == 0):
            to_check = []
        # if there is 1 word we use it 
        elif (len(useful_words) == 1):
            # the records to compare to all other records containing it
            list1 = set(title_index[useful_words[0][0]]['docs']) #list(set(title_index[useful_words[0][0]]['docs']))
            to_check = [x for x in list1 if coll_keys.index(x) > i]
            del(list1)
        else:
            # if there are 2 words or more, we identify the two words closest to median df
            below_median = (len(useful_words)-1)//2
            # the records to compare must contain at least one of these words
            set1=set(title_index[useful_words[below_median][0]]['docs'])
            set2=set(title_index[useful_words[below_median+1][0]]['docs'])
            orlist = set1 | set2 #list(set1 | set2)
            del(set1)
            del(set2)
            to_check = [x for x in orlist if coll_keys.index(x) > i]
            del(orlist)
    else:
        if (len(useful_words) < 2):
            to_check_range = range( i+1, len(coll_keys) )
            to_check = [coll_keys[i] for i in to_check_range]
        else:
            # if there are 2 words or more, we identify the words closest to the median df
            below_median = (len(useful_words)-1)//2
            # the records to compare must contain at least one of these words
            set1=set(title_index[useful_words[below_median][0]]['docs'])
            set2=set(title_index[useful_words[below_median+1][0]]['docs'])
            orlist = set1 | set2 #list(set1 | set2)
            del(set1)
            del(set2)
            to_check = [x for x in orlist if coll_keys.index(x) > i]
            del(orlist)
        
    for kk in to_check:
        recid_pairs.append( (coll_keys[i],kk) )
    del(to_check)

    return recid_pairs

###############
### Logging ###
###############

def write_message(log, message):
    with open(log,'a') as f:
        f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S -> ") + str(message) +'\n')
