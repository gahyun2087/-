# PCR Primer Analyzer

Python으로 구현한 PCR 프라이머 분석 프로그램 (기능 1: 기존 프라이머 분석)

## 기능

- 표적 DNA, Forward/Reverse primer 입력 및 서열 유효성 검사 (A/T/G/C 외 문자 검출)
- 프라이머 길이, GC 함량, Tm 계산
  - Tm은 단순 계산식(Wallace rule)과 Primer3-py의 열역학적 계산(Nearest-Neighbor)을 함께 제공
  - 실제 분석에는 Primer3-py의 열역학적 계산 값을 기본으로 사용
- 표적 DNA 내 프라이머의 **모든** 결합 위치 탐색 (reverse primer는 역상보서열로 변환하여 탐색)
- 결합 위치 조합을 검토하여 **올바른 방향(Forward가 앞, Reverse가 뒤에서 마주보는 배치)** 인 경우에만 PCR 산물 길이 계산
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

## 알려진 한계

- Primer3-py의 `calc_heterodimer`/`calc_homodimer`는 서열 전체를 기준으로 가장 안정적인 결합 구조 하나만 계산한다. "3' 말단 부위만을 따로 떼어 다이머 여부를 판정하는 기능"은 Primer3-py가 별도로 제공하지 않으므로 본 프로그램에도 구현하지 않았다.
- 결합 위치 탐색은 완전 일치(exact match)만 지원하며, 미스매치를 허용한 결합은 탐색하지 않는다.

## 진행 상황

- [x] 서열 입력 및 유효성 검사
- [x] 프라이머 길이 / GC 함량 / Tm 계산
- [x] 모든 결합 위치 탐색 및 방향 검증을 포함한 PCR 산물 길이 계산
- [x] Forward/Reverse Tm 차이(ΔTm) 계산
- [x] Hairpin / Self-dimer / Heterodimer 구조 분석
- [x] 결과 요약 및 기준 충족 여부 출력 ([양호]/[주의] 표시)
- [ ] 기능 2: Primer3-py 기반 신규 프라이머 후보 생성 (예정)
