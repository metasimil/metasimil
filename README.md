# MetaSimil

MetaSimil is a metadata similarity analysis framework written in pure Python. 

- Author: Jan Krause
- Current Version: 0.1 (alpha), as the framework is not fully tested the use of MarcXimiL is recommended for production.
- First release: 2016-08-28
- Licence: [GNU GLPv3](http://www.gnu.org/copyleft/gpl.html)

## Features

- Duplicates detection or similar records suggestion,
- Single or multiple collections analysis, 
- Input format: [DublinCore](http://dublincore.org/schemas/xmls/), [MarcXML](http://www.loc.gov/standards/marcxml/), [CSV](https://en.wikipedia.org/wiki/Comma-separated_values), or Python objects (when used as a library),
- Live analysis within a data repository,
- Based on efficient information retrieval algorithms, 
- Can handles hundred of thousands of records on a laptop,
- Multi-processing.

# Tutorial an examples

See this [general example](./Examples/AnalysisScripts/AnalyseOneCollection.py).

## Origin and credits

MetaSimil is an evolution of [MarcXimiL](http://marcximil.sourceforge.net/), written by Alain Borel and Jan Krause and which is also licensed under GPLv3. Notable changes:

- extension of input formats: MarcXML, DublinCore, CSV and Pyhon objects,
- easily usable in Python scripts (just import the library),
- selection of the very best and most useful algorithms of MarcXimiL,
- improvement of some algorithms (notably date and authors similarity computation)
- improvement in speed and memory consumption.

