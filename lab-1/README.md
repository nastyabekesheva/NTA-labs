# NTA lab 1
Script to factorize numbers. Uses trial divison, rho-Pollard and Brillhart-Morrison algoritms and Solovay-Shtrassen test.

## Report

You can read my report in `/report` folder. It has time stats for factorization of a few numbers and description of problems I encountered while doing this lab.

## Usage
pull instance from dockerhub 

```
docker pull bekeshevaaa/nta-lab-1:0.6
```

run specifing parameters [number to factorize]

```
docker run bekeshevaaa/nta-lab-1:0.6 python3 script.py num
```
