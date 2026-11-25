# PCR Primer Analyzer

Python으로 구현한 PCR 프라이머 분석 프로그램 (기능 1: 기존 프라이머 분석)

## 기능

- 표적 DNA, Forward/Reverse primer 입력 및 서열 유효성 검사
  - 공백, 탭, 줄바꿈 등 모든 공백 문자 제거 후 대문자로 통일
  - A/T/G/C 외 문자, 빈 서열 입력 시 오류 메시지 출력
- 프라이머 길이, GC 함량, Tm 계산
  - Tm은 단순 계산식(Wallace rule)과 Primer3-py의 열역학적 계산(Nearest-Neighbor)을 함께 제공
  - 실제 분석에는 Primer3-py의 열역학적 계산 값을 기본으로 사용
  - Tm 계산에 사용한 농도 조건(Na+, Mg2+, dNTP, 프라이머 농도)을 결과에 명시 (현재는 고정 기본값 사용)
- 표적 DNA 내 프라이머의 **모든** 결합 위치 탐색 (reverse primer는 역상보서열로 변환하여 탐색)
- PCR 산물 길이 계산 시, **Forward 결합 구간이 끝난 뒤 Reverse 결합 구간이 시작하는 경우**(겹치지 않고 서로 마주보는 올바른 배치)만 유효한 조합으로 인정
  - 조건을 만족하는 **모든** 조합을 번호와 함께 나열 (가장 짧은 산물을 자동으로 추천하지 않음)
  - 유효한 조합이 없으면 겹침/방향 문제 등 이유를 안내
- Forward/Reverse primer의 **Tm 차이(ΔTm)** 계산 및 2℃ 이내 여부 확인
- Hairpin, Self-dimer, Heterodimer 구조 분석 (ΔG, Tm)
- 일반적인 프라이머 설계 가이드라인과 비교한 결과 요약 출력 ([양호]/[주의] 표시)

## 설치

```bash
pip install -r requirements.txt
```

## 실행

```bash
python primer_analysis.py
```

실행 후 안내에 따라 표적 DNA, Forward primer, Reverse primer 서열을 순서대로 입력하면 된다.

## 사용 라이브러리

- [Biopython](https://biopython.org/): DNA 서열 처리 및 역상보서열 계산
- [primer3-py](https://github.com/libnano/primer3-py): Tm 계산 및 Hairpin/Self-dimer/Heterodimer 등 열역학적 구조 분석

## 참고 기준 (일반적인 설계 가이드라인)

| 항목 | 권장 범위 |
|---|---|
| 프라이머 길이 | 18~24 bp |
| GC 함량 | 40~60 % |
| Tm (열역학적 계산) | 50~65 ℃ |
| Forward/Reverse Tm 차이(ΔTm) | 2 ℃ 이내 |
| 구조 ΔG (Hairpin/Self-dimer/Heterodimer) | -9 kcal/mol 초과 권장 (참고용 경고 기준, 절대 기준 아님) |

위 수치는 일반적으로 통용되는 참고 기준이며, 절대적인 합격/불합격 기준이 아니다. 본 프로그램의 출력 결과는 실제 PCR 실험의 성공을 보장하지 않는다.

## Tm 계산 조건 (기본값)

현재 프로그램은 아래 고정값으로 Tm을 계산하며, 실제 PCR 반응 조건을 입력받는 기능은 아직 없다. 계산 결과 하단에 이 조건이 항상 함께 표시된다.

| 항목 | 기본값 |
|---|---|
| 1가 양이온(Na+/K+ 등) 농도 | 50 mM |
| Mg2+ (2가 양이온) 농도 | 0 mM |
| dNTP 농도 | 0 mM |
| 프라이머(DNA) 농도 | 50 nM |

## 알려진 한계

- Primer3-py의 `calc_heterodimer`/`calc_homodimer`는 서열 전체를 기준으로 가장 안정적인 결합 구조 하나만 계산한다. "3' 말단 부위만을 따로 떼어 다이머 여부를 판정하는 기능"은 Primer3-py가 별도로 제공하지 않으므로 본 프로그램에도 구현하지 않았다.
- 결합 위치 탐색은 완전 일치(exact match)만 지원하며, 미스매치를 허용한 결합은 탐색하지 않는다.
- FASTA 형식 입력은 아직 지원하지 않는다 (일반 DNA 서열 문자열만 입력 가능).
- 실제 PCR 반응 조건(농도 등)을 직접 입력받는 기능은 아직 없으며, 고정된 기본값으로만 Tm을 계산한다.

## 진행 상황

- [x] 서열 입력 및 유효성 검사 (모든 공백 문자 제거 포함)
- [x] 프라이머 길이 / GC 함량 / Tm 계산 (계산 조건 명시)
- [x] 모든 결합 위치 탐색 및 겹침/방향 검증을 포함한 PCR 산물 길이 계산 (모든 유효 조합 표시)
- [x] Forward/Reverse Tm 차이(ΔTm) 계산
- [x] Hairpin / Self-dimer / Heterodimer 구조 분석
- [x] 결과 요약 및 기준 충족 여부 출력 ([양호]/[주의] 표시)
- [ ] 기능 2: Primer3-py 기반 신규 프라이머 후보 생성 (예정)
