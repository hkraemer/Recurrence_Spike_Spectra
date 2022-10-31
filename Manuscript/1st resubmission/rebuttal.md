Rebuttal letter
===============

We thank the reviewers for their careful reading of our manuscript and their useful comments. We have considered all of them and improved the manuscript accordingly. Our responses to all reviewer comments (listed below in italic) in the following. 

Reviewer 1
----------

_Introduction: a clear formulation of the mathematical problem, the state-of-the-art in terms of the existing literature, a summary of the solution proposed and of its pros and cons are totally lacking._

We thank for this comment, but do not agree completely with this criticism. In our opinion, the problem is well stated and the outline of the manuscript is summarised. But the literature review can be improved, which with the help of the suggestions of the reviewer 2 we finally did. "a summary of the solution proposed" as well as "pros and cons" do not belong to the introduction. These are provided in the conclusion of the manuscript. We added some further introductory sentences to the beginning of the introduction.

_Equation (1) is poorly introduced. Although the authors cite a reference ([1]), the x_i, x_j vectors, as well as the parameters d and N, are not explained. The "trajectory x_i", for example, is introduced only at the end of page 1. Moreover, R_{i,j}(\epsilon) is not explicitly stated and explained._

Indeed, out intention was to have a short but clear introduction into the recurrence plot method, which reviewer 2 obviously liked very much. Nevertheless, we have now added an explanation on the recurrence plot definition. An unfamiliar reader should now better understand the method.

_(lines 20-21) "Furthermore, the correlation structures of higher dimensional spaces can be resolved in the recurrence-derived Fourier-spectrum." What do the authors mean with "correlation structures of higher dimensional spaces"? And what do they mean with the term "resolve"?_

TBD

_Figs. 2C, 2F: the "Inter Spike Spectrum of perfect DC" is equal to the "Inter Spike Spectrum of randomized DC", although the time domain signals (Figs. 2A, 2D, respectively) are not. This might prove that the new tool introduced in the manuscript is not invertible: Definitely a major issue._

TBD

_(beginning of Sec. 2) "Let s(t_i) be the normalized signal we want to transform...". It is unclear how s(t_i) is obtained (i.e. "normalized") from x(t)._

TBD

_The way in which basis functions are defined is confusing._

TBD

_The authors state that they "either use" the LASSO or the STLS methods to obtain β. It is not clear when the authors prefer one method over the other one, and why. Moreover, although a direct comparison between the methods, namely the different behavior of the regularization parameter, is shown in Fig. A3, a thorough comparison between the two methods is lacking: Figs. 2, 4 and 5 appear to be obtained by using only the LASSO method. One of the few mentions of the STLS method is in the Discussion section, where it is stated that (lines 236-237) "[...] the two different sparse regression algorithms [...] yield different results for the same desired ρ". How do the authors support this statement? There is no clear way for a reader to either compare the two methods, or to understand why the authors did not consider a single method from the very beginning._

TBD

_At line 51 the authors introduce the term "loading": The authors should explain the meaning of this term._

TBD

_With regard to the analysis of power grid frequency time series (Sec. 3.3, Fig. 6) The criterion followed by the authors to evaluate the importance of the spectral peaks is unclear._

TBD

_Overall, the details of the power spectrum estimation procedure are not fully described: in particular, no information on windowing is provided. It is thus hard to evaluate the significance of the "peak splitting" discussed by the authors in Sec. 3.3 and concerning Fig. 6(A): could this splitting be due to spurious effects such as spectral leakage? It is worth noting that the information on windowing is instead provided in the case of the "evolutionary" spectra of Earth's orbit data._

TBD

_The claimed robustness to noise is evaluated by showing the results of the method applied to the Roessler system with 5% additive Gaussian white noise (Fig. A2). While the results indeed appear to be unchanged with respect to the noiseless case (Fig. 4), a single example is hardly sufficient to claim that the method is "robust to noise". At which signal-to-noise ratio does the method fail to provide reliable results?_

TBD

_The authors unduly use the term "powerspectrum". However, this word does not exist in English. They might consider, for example, power spectral density._

We do not agree with this criticism. "Powerspectrum" (or "power spectrum") is common term in data analysis, see, e.g., <https://mathworld.wolfram.com/PowerSpectrum.html>. Nevertheless, we have modified the term to "power spectral density".

_The abbreviation "DC" for "Dirac comb" is introduced four times (lines 27, caption of Fig. 2, beginning of Sec. 2, line 248)._

TBD

_The abbreviation "s.t." in Eq.(6) is unclear. Does it stand for "such that"? Why not using "where"?_

"s.t." means "such that" or "so that" and is frequently used in equations. Nevertheless, we have now explicetely given "such that" (which is correct, and not "where", because the following part is the condition for variable $n$).

_(lines 241-244) "This is not a drawback of the proposed method, but rather a drawback of the particular application method that we have heavily used in this article and which was the main motivation for developing the proposed method." This sentence is unclear: which "application method" are the authors referring to?_

TBD

 _Some typos:
 in page 2, Wiener-Khinchim ---> Wiener-Khinchin;_
 
 We corrected this typo.
 
 _in the Conclusion, "Wee chose LASSO..." (unless Wee is a surname)._

 We corrected this typo.
 

Reviewer 2
----------

_Starting a detailed assessment from the beginning, it is necessary to emphasize the comprehensive and impeccable introduction to the issue of the construction of a novel kind of power spectrum, which is the inter-spike spectrum._

We have extended the introduction and hope it is now more clear.

_Moreover, the Authors explained the algorithm in great detail (in numerous pictures). By the way, also with outer -spike. Literature examples below review._

We have added some further references on spike-train analysis and, in particular, on spike-train power spectra. Thanks for the literature suggestions, which were really helpful!

_The illustrations of tau-recurrence rate based spectrum (Fig. 1) and the transformation diagrams of the series of Dirac delta function (Fig. 2), including the two proposed inter-spikes spectrum, deserve emphasis. A small remark - regarding verse 35: what will be with stochastic resonance, because we are dealing with impulse disturbances.?_

TBD

_Chapter 2 presents the signal decomposition into set of appropriate basis functions using Dirac comb. Here, too, we can distinguish a graphic representation of the decomposition procedure (Fig. 3). Nevertheless, the Authors could cite a position in the literature that discusses other methods of decomposition._

We took the reviewer's advice and think this is a good point. Hence, we extended the Methods-section and added the corresponding references.

_A small note on Appendix C: It would be interesting to present an inter-spike spectra also for one of the nonlinear systems from Duffing, Mackey-Glass or Chen systems._

We thank the reviewer for this suggestion. In fact, we were computing numerous inter-spike spectra for different systems. But regarding the manuscript, we decided to present a variety of possible applications, which includes the analysis of exemplary nonlinear system (like the Roessler we showed or systems you suggest) as well as "real-world" data (power grid frequency data) and data stemming from a geophysical context (orbital data). Each of these applications could, of course, be expanded or other applications could be shown, but here we must concentrate on presenting the method itself as concisely as possible while keeping the size of the manuscript within tolerable limits. Moreover, when analyzing a nonlinear system we would potentially be interested in different dynamical setups, which would, in turn, lead to a variety of aspects we would additionally need to consider and cover with the proposed method. The idea of this paper is to outline the proposed method and provide the code/package for its use. We would certainly be pleased if researchers decided to publish a thorough study of a number of different nonlinear systems using this method.