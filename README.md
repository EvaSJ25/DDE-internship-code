# DDE-internship-code

**Table of contents**
1. [Introduction](#introduction)
2. [Barycentric interpolation](#barycentric-interpolation)
3. [DDETools](#DDETools)
4. [Stability using Breda et al. 2009 method](#stability-using-Breda-et-al.-2009-method)
5. [Mackey-Glass equation demo](#mackey-glass-equation-demo)
6. [Neuron demo](#neuron-demo)
7. [Single-delayed inverted pendulum demo](#single-delayed-inverted-pendulum-demo)
8. [References](#references)

## Introduction

This repository contains demos and code for some of the tools required for delay differential equation (DDE) analysis. These functions include finding equilibria (or steady-states (stst)); finding stability (using two different methods); finding Hopf and fold bifurcations and also the continuation of these bifurcation over a 2-parameter space. To demonstrate how my code can be implemented I have included demos in Jupyter notebooks for three examples: Mackey-Glass equation, a single-delayed inverted pendulum problem and an example of modelling interaction between two neurons (the system of this is the same as the neuron example given for MATLAB's DDE-Biftool).

## Barycentric interpolation

One method of stability finding, the Breda et al method (Breda et al. 2009), is built on barycentric interpolation. The Jupyter notebook contained in this folder talks through the background information on barycentric interpolation before going through some examples of how to use the functions I created for use in barycentric interpolation.

## DDETools

This folder contains the file for the module DDETools - a module I created that contains the functions I developed for DDE analysis in Julia. The 'shared' folder contains the developed functions

## Stability using Breda et al. 2009 method

This folder contains a Jupyter notebook that discusses the theory behind the method put forward by Breda et al. for finding DDE stability in the paper (Breda et al. 2009). It introduces and explains the functions I have created to implement this method before it walks through examples to aid user understanding.

## Mackey-Glass equation demo

This folder contains the functions that define the Mackey-Glass equation (and its delay) and are in a form where they can be used in my DDE analysis functions. It also contains a Jupyter notebook that introduces the reader to the Mackey-Glass equation and then shows how the developed functions in DDETools can be used to analyse this DDE system.

## Neuron demo

This folder contains the functions that define a system to model the reaction between two neurons (and its delays). This system is the same as the neuron demo for MATLAB's DDE-Biftool (*DDE-BIFTOOL demo 1 - Neuron*, 2014). The folder also contains a Jupyter notebook that introduces the reader to the neuron DDE system and then shows how the developed functions in DDETools can be used to analyse this DDE system.

## Single-delayed inverted pendulum demo

This folder contains the functions that define a simplified system for a single-delayed inverted pendulum problem (and its delays). The folder also contains a Jupyter notebook that introduces the reader to the neuron DDE system and then shows how the developed functions in DDETools can be used to analyse this DDE system.

## References

1. *DDE-BIFTOOL demo 1 - Neuron* (12/2014). Available at https://ddebiftool.sourceforge.net/demos/neuron/html/demo1.html (Accessed: 29 August 2025)
2. Breda D., Maset S., Vermiglio R. (2009). 'TRACE-DDE: a Tool for Robust Analysis and Characteristic Equations for Delay Differential Equations', *Lecture Notes in Control And Information Sciences*, volume 388, pp 145-155. Available at: http://dx.doi.org/10.1007/978-3-642-02897-7_13