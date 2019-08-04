# Splice_Deep
Splice2Deep: Deep Learning Model for Recognition of Splice Sites in Genomic DNA.

The accurate identification of the exon and intron boundaries is critical for the correct annotation of genes with multiple exons contained in a primary genomic DNA sequence. These boundaries are demarcated by the donors and acceptors splice sites (SS). Deriving accurate computational models to predict the SS are useful for functional annotation of genes and genomes, and for finding alternative SS that are frequently associated with different diseases. Although different models have been proposed for the computational prediction of SS, there is a need to improve the accuracy further to enable a reliable annotation. Moreover, models are often derived and tested by using the same organism, therefore essentially not providing information about broader generalization. 
We developed the SpliceDeep model, an ensemble of deep convolutional neural networks, to automatically extract the features from the DNA and being capable of generalizing results to other organisms for which it was not trained. The performance of the model is evaluated on the recognition of SS of different organisms: Homo sapiens, Oryza japonica, Arabidopsis thaliana, C.elegans, and Drosophila melanogaster. In addition to improved accuracy compared to the state-of-the-art tools, the SpliceDeep results obtained from the cross-organism validation demonstrated that the model correctly identified conserved genomic elements that may allow users to annotate SS in new genomes by choosing the taxonomically most close model.

### Requirements
  - Models run on linux machine.
  - Anaconda Python 2.7 or later.
  - keras.
    
### Usage
We propose five models to recognize SS in different organisms: Homo sapiens 'hs', Oryza japonica 'oriza', Arabidopsis thaliana 'at', C.elegans 'c_elegans', and Drosophila melanogaster 'd_mel'. Our models take DNA seqences with 602 length (SS in 300 and 301 positions) as input, performs feature extraction and feature selection using DL from the flanking regons, and predict whether the given sequence represents a true/false SS using an artificial neural network (NN) binary classifier.    

To run an Acceptor model:

```
  $ python Splice_Deep_Acceptor.py org='pass a shortcut of organism of intrest: hs,at,c_elegans,d_mel,oriza' fname= 'pass fasta file'
```
  # e.g. to use a Homo sapiens model: (Note: AcSS_test is provided in Data folder and contain 2000 positive/negative seqences)
```  
  $ python Splice_Deep_Acceptor.py org='hs' fname='AcSS_test.fa' 
    
```
Predections will be stored in 'splicedeep_AcSS_output.txt' file. 


To run a Donor model:

```
  $ python Splice_Deep_Donor.py 'pass a shortcut of organism of intrest: hs,at,c_elegans,d_mel,oriza' fname= 'pass fasta file'
  ```
  # e.g. to use a Homo sapiens model: (Note: DoSS_test is provided in Data folder and contain 2000 positive/negative seqences)
```  
  $ python Splice_Deep_Donor.py org='hs' fname='DoSS_test.fa'
  
```
Predections will be stored in 'splicedeep_DoSS_output.txt' file.  
  

