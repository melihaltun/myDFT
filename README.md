# myDFT
Real and complex Discrete Fourier Transform and inverse Discrete Fourier Transform functions

Even though FFT is faster than DFT, sometimes memory is the limiting factor rather than time. Moreover, for small arbitrary size inputs, the computational overhead of buffer allocations may negate the speed advantage of FFT. DFT also has the advantage of being able to calculate specific frequency bins without calcuating the whole spectrum. In that particular use case DFT computational complexity is O(N) which beats the O(NlogN) complexity of FFT. In such cases DFT might be better suited to an application than FFT.

This code has several functions that can perform real and complex DFT/iDFT operations. There also functions that target specific frequency bins and only calculate results for them. Before the operations a DFT object must be declared and its buffers must be allocated with the set_dft_instance() function. After the operations delete_dft_instance() function is required to release the allocated memory.  
