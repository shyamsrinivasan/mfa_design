\ Model mfa_design
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  1.224744871391589 mut_select_cont[0]
   + 1.224744871391589 mut_select_cont[1]
   + 2.449489742783178 mut_select_cont[2]
Subject To
 c1[0,0,1]: - mut_select[0] + mut_select_cont[0] <= 0
 c1[1,0,2]: - mut_select[0] + mut_select_cont[1] <= 0
 c1[2,1,2]: - mut_select[1] + mut_select_cont[2] <= 0
 c2[0,0,1]: - mut_select[1] + mut_select_cont[0] <= 0
 c2[1,0,2]: - mut_select[2] + mut_select_cont[1] <= 0
 c2[2,1,2]: - mut_select[2] + mut_select_cont[2] <= 0
 c3[0,0,1]: mut_select[0] + mut_select[1] - mut_select_cont[0] <= 1
 c3[1,0,2]: mut_select[0] + mut_select[2] - mut_select_cont[1] <= 1
 c3[2,1,2]: mut_select[1] + mut_select[2] - mut_select_cont[2] <= 1
 c4: mut_select[0] + mut_select[1] + mut_select[2] <= 5
Bounds
 mut_select_cont[0] <= 1
 mut_select_cont[1] <= 1
 mut_select_cont[2] <= 1
Binaries
 mut_select[0] mut_select[1] mut_select[2]
End
