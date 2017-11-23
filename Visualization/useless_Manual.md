#16s following analysis script description

##1. Rarefaction curve
```
Script name: draw_PD.py
requiremnts: 
python2.7
scikit-bio == 0.4.2
plotly == 2.0.12
pandas
```

example:
python draw_PD.py -t /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otus.tree -i /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/analysis/otu_raw_filterd_e5.txt -m /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/sample_info.csv -o /home/liaoth/data2/hospital_seven/PG17H11010002Amp01/drawn/


##2. MDS plus MST & KNN
```
Script name: DRAW_MDS_plus.py
requiremnts: 
python2.7
scikit-learn == 0.18.1
networkx == 1.11
plotly == 2.0.12
pandas
```

Majorly, it just need to modify comment #Config below code.

If you doesn't need to rename or group display dots. You can just assign `None` to `metadata,col_need,otherinfo_col`

Output dir you doesn't need to build it.
##3. Stack Bar for species distribution
```
Script name: draw_stack_bar.py
requiremnts: 
python2.7
scikit-bio == 0.4.2
plotly == 2.0.12
pandas
```
Majorly, it just need to modify comment #Config below code.


##4. Bortuta to extract distinct feature
```
Script name: construc_bortuta_input.py | draw_bortuta_result.py
requiremnts: 
python2.7
scikit-bio == 0.4.2
plotly
pandas
```

lazy to write. short and i believe you can understand. Inside have two plot are draw.