# Import Libraries

import numpy as np
import pandas as pd
import streamlit as st
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import subprocess
import os
import base64

# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# Lipinski Descriptor calculator
def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData= np.arange(1,1)
    i=0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1

    columnNames=["MW","NumHDonors","NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)

    return descriptors

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('aromtase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Page title
st.markdown("""
# Aromtase Bioactivity Prediction App
This app allows you to predict how effective a compound is at inhibting the `Aromatase` enzyme. `Aromatase` is a drug target for Breast Cancer.
**Credits**
- App built in `Python` + `Streamlit` by [Gowrav Mannem](https://www.linkedin.com/in/gowrav-mannem-830896218/).
- Adopted from Chanin Nantasenamat's(AKA [Dataprofessor](https://github.com/dataprofessor)) [Youtube Tutorial](https://www.youtube.com/watch?v=jBlTQjcKuaY&t=4960s)
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/).
- Lipinski Descriptors found using [RDKit](https://www.rdkit.org/).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")
# User starts prediction process
if st.sidebar.button('Predict'):
    # reading input data
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    # printing out input data
    st.header('**Original input data**')
    st.write(load_data)

    # animation
    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated Fingerprint descriptors**')
    fp_desc = pd.read_csv('descriptors_output.csv')
    st.write(fp_desc)
    st.write(fp_desc.shape)

    # Read in Lipinski descriptors
    lipinski_df=lipinski(load_data)
    st.header('**Calculated Lipinski Descriptors**')
    st.write(lipinski_df)
    st.write(lipinski_df.shape)

    # combining the two dataframes
    aromatase_XY = pd.concat([fp_desc, lipinski_df], axis=1).reindex(fp_desc.index)
    # Read descriptor list used in previously built model
    st.header('**Combining the two dataframes and dropping Low variance features**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = aromatase_XY[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')
