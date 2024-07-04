# Part assembly 자동화를 위한 설계

- 부품/모듈 조합을 설계하고 genbank 및 fasta 파일 생성
- 설계는 엑셀로 수행하고 해당 파일을 입력으로 genbank 및 fasta 파일을 생성해주는 파이썬 코드 활용
- 설계는 다음 네 단계로 이루어짐
  - 부품 정보 관리 (PartDB-kribb.xlsx)
  - 부품 단위 Golden gate assembly 호환 설계 (part_preparation.xlsx)
  - 부품 어셈블리 설계 (part_assembly.xlsx)
  - 모듈 어셈블리 설계 (module_assembly.xlsx)

## Part DB
-   PartDB-kribb.xlsx
-   PlateA_xx, PlateB_xx 유지관리 중지
-   향후 stock 정보용 별도 DB 유지관리 필요 

## 부품 입출고 자동화를 위한 프로세스 정립

-   베타 장비 중 Incubator 저장 및 관련 부품 자동 입출고 프로세스 고민 필요
-   임의의 DNA 부품 (일정 농도유지, 정기적으로 관리)
-   screw tube 저장
-   각 위치별 mapping 정보 저장고 sw 같이 활용 필요


## Part assembly를 통한 모듈 준비 과정 (MVA 예시)

- Variance region (출처: yeast vegas assembly)
  - VA16,TATCGCGGGTGCGTGCATCGACAAGCCATGCCCACCTTCTGGTCGATTGGGCTGGCG,NA,spacer,NA
  - VA3,GGAGGTACTGGCCTAGCGTCGTGGCCCGGGAGAGACAGTTTAGTAGTGACTCGCGGC,NA,spacer,NA
  - VA4,TTGGCGTTAATTGTAGCTTATTTCCCGCCCTGTGATTGAGGCGGGATGGTGTCCCCA,NA,spacer,NA
  - VA5,GACTAAGACTCTGGTCACGGTTCAGAAGTGGACGATGCATGTCGTCGGGCTGATAGA,NA,spacer,NA
  - VA6,TGCACGGCGCTAGGTGTGATATCGTACACTTGGGAGAAGTCAGATACGATTGCGGCT,NA,spacer,NA
  - VA7,TAGCGGCGCCGGGAAATCCAGCATATTCTCGCGGCCCTGAGCAGTAGGTGTCTCGGG,NA,spacer,NA
  - VA8,GAGTCTACGTTACACCTGAACTCGCATGTCTGGGGTTGTGGTCAGGCCTTGTCAATT,NA,spacer,NA
  - VA9,GCGTACTGGCCGCCCGGGCCTGATGTGGCCGTCCTATTAGCATTGTACACCCTCATT,NA,spacer,NA
  - VA2,TGACGCTTGGATGCGTGACCCCGTACGTCATGACCCGTCATGGGTATGTAAGCGAAG,NA,spacer,NA

- Overhang 
  - O1,GCCT,NA,overhang,NA vector - promoter 
  - O2,CTTT,NA,overhang,NA promoter - rbs
  - O3,GCAG,NA,overhang,NA rbs - cds
  - O4,CTAA,NA,overhang,NA cds - term
  - O5,TCAC,NA,overhang,NA term -vector

- 모듈 예시: VA16 - Mvak1 - VA3  construct는 아래와 같음

      (VA16 duplex 준비)
      5' - VA16 - O1 - BsaI - 3' 
      5' - BsaI - O1 - promoter - O2 - BsaI - 3'
                      5' - BsaI - O2 - rbs - O3 - BsaI - 3'
                                5' - BsaI - O3 - mvak1 - O4 - BsaI - 3'
                                              5' - BsaI - O4 - term - O5 - BsaI - 3'
                                                          5' - BsaI - O5 - VA3 - 3'

- 순차적으로 다음 모듈들 제작 (genbank 파일 제작)
      VA16 - MvaK1 - VA3 
      VA3 - MvaD - VA4
      VA4 - MvaK2 - VA5
      VA5 - Idi - VA6
      VA6 - MvaE(bsaI del) - VA7
      VA7 - MvaS - VA8
      VA8 - IspA - VA2 

## 모듈 어셈블리를 통한 모듈 조합 준비

            VA16 - mvak1 - VA3 - oh - BsaI
BsaI - oh - VA3 - MvaD - VA4 - oh - BsaI
BsaI - oh - VA4 - MvaK2 - VA5 

            VA5 - Idi - VA6 - oh - BsaI
BsaI - oh - VA6 - MvaE(bsaI del) - VA7 - oh - BsaI
BsaI - oh - VA7 - MvaS - VA8 - oh - BsaI
BsaI - oh - VA8 - IspA - VA2

Golden Gate Assembly 

VA16 - VA3 - VA4 - VA5
VA5 - VA6 - VA7 - VA8 - VA2

Vector 

VA16 - lycopen operon - Vector - VA2

VA16 --- VA3 --- VA4 --- VA5 

PCR primer? 






