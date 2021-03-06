#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 05-Jun-2012
#
# Build for R2.13.2, 32bits
#-------------------------------------------------------------------------------

#- Scenario descriptions
Fmsy          <- 0.25


LTMP          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.2e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP")                                                                                                                                                                                                                                                                                                                                                                                                
                                )
                                
LTMPAlpha          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.1,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.11,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.12,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.13,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.14,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.15,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.16,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.17,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.18,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.19,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.20,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=1,
                                 alpha=0.21,
                                 stabilityBreak="Btrigger",scen="LTMPAlpha")                                                                                                                                                                                                                                                                                                                                                                                                
)

                                
                     
LTMP2mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP2mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP2mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )


                     
LTMP2.4mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.4e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP2.4mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )                                              
                                
                                
                     
LTMP2.6mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.6e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP2.6mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )                                 
                                

                     
LTMP2.8mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=2.8e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP2.8mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )                                                                 
                                
                                
                                
                     
LTMP3.0mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.0e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP3.0mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )                                  
                                
                                
LTMP3.2mt          <- list(opt1=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.01,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt2=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.02,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt3=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.03,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt4=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.04,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt5=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.05,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt6=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.06,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt7=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.07,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt8=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.08,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt9=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.09,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt10=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.1,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt11=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.11,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt12=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.12,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt13=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.13,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt14=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.14,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt15=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.15,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt16=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.16,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt17=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.17,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt18=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.18,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt19=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.19,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt20=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.20,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt"),
                      opt21=list(Blim=1.84e6,Bpa=2.36e6,Btrigger=3.2e6,
                                alpha=0.21,
                                stabilityBreak="Btrigger",scen="LTMP3.2mt")                                                                                                                                                                                                                                                                                                                                                                                                
                                )                                     