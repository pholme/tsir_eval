# Evaluating the fast event-driven algorithm for the SIR model on temporal networks

This code compares the even-driven implementation of the SIR model on temporal networks available at https://github.com/pholme/tsir with a straightforward implementation. It generates temporal networks—binomial random graphs with n nodes, z expected average degree, c contacts per edge on average.

To compile the code, make a directory for the object files `mkdir o` and compile the code `make`.
Given network parameters (n, z, c) and SIR parameters. You can run it as:

```python3 tsir_compare.py [n] [z] [c] [beta] [nu]```

The output gives the relative speed-up (how many times faster the event-driven model is to the straightforward implementation), the p-value of a Mann-Whithey U-test of the output and the average outbreak size.
