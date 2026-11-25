"""
primer_analysis.py

목적: 사용자가 입력한 표적 DNA 서열과 Forward/Reverse primer를 분석하는 프로그램
- 서열 유효성 검사
- 프라이머 길이 / GC 함량 / Tm 계산
- 표적 DNA 내 모든 결합 위치 탐색 및 올바른 방향의 PCR 산물 길이 계산
- Forward/Reverse Tm 차이(ΔTm) 확인
- Hairpin / Self-dimer / Heterodimer 구조 분석
- 일반적인 설계 기준과 비교한 결과 요약 출력

사용 라이브러리: Biopython, primer3-py
주의: 본 프로그램의 판정은 일반적으로 통용되는 설계 가이드라인과의 비교일 뿐이며,
      실제 PCR 실험의 성공 여부를 보장하지 않는다.

참고: Primer3-py의 calc_heterodimer / calc_homodimer는 두 서열 전체를 기준으로
      가장 안정적인 결합 구조 하나를 계산한다. "3' 말단 부위만을 따로 떼어
      다이머 여부를 판정하는 기능"은 Primer3-py가 별도로 제공하지 않으므로
      이 프로그램에서도 구현하지 않는다. (검증되지 않은 기능을 구현한 것처럼
      표시하지 않기 위함)
"""

import primer3
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# 1. 입력 및 서열 유효성 검사
# ---------------------------------------------------------------------------

def clean_sequence(seq):
    """서열을 받아서 대문자로 바꾸고, 공백/줄바꿈을 제거한다."""
    seq = seq.strip()
    seq = seq.upper()
    seq = seq.replace(" ", "")
    return seq


def is_valid_dna(seq):
    """서열이 A, T, G, C 로만 이루어져 있는지 검사한다."""
    allowed_bases = set("ATGC")
    return set(seq).issubset(allowed_bases)


def get_sequence_input(prompt_text):
    """사용자로부터 서열을 입력받고, 정리 및 검증까지 수행한다."""
    raw_seq = input(prompt_text)
    seq = clean_sequence(raw_seq)

    if seq == "":
        print("오류: 서열이 입력되지 않았습니다.")
        return None

    if not is_valid_dna(seq):
        print(f"오류: '{raw_seq}' 에 A, T, G, C 이외의 문자가 포함되어 있습니다.")
        return None

    return seq


# ---------------------------------------------------------------------------
# 2. 기본 특성 계산 (길이, GC 함량, Tm)
# ---------------------------------------------------------------------------

def calc_length(seq):
    """서열 길이를 반환한다."""
    return len(seq)


def calc_gc_content(seq):
    """GC 함량(%)을 계산한다."""
    g_count = seq.count("G")
    c_count = seq.count("C")
    gc_percent = (g_count + c_count) / len(seq) * 100
    return round(gc_percent, 2)


def calc_tm_simple(seq):
    """
    단순 계산식(Wallace rule)으로 Tm을 추정한다.
    - 짧은 프라이머(대략 14bp 이하)에 적합한 근사식
    - 반응 조건(염 농도 등)을 전혀 반영하지 않음
    """
    a_count = seq.count("A")
    t_count = seq.count("T")
    g_count = seq.count("G")
    c_count = seq.count("C")
    return 2 * (a_count + t_count) + 4 * (g_count + c_count)


def calc_tm_thermo(seq, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=50.0):
    """
    Primer3-py의 열역학적 계산(Nearest-Neighbor 방식)으로 Tm을 계산한다.
    실제 반응 조건(Mg2+ 농도 등)을 파라미터로 반영할 수 있어 Wallace rule보다 정확도가 높다.
    """
    tm = primer3.bindings.calc_tm(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
    )
    return round(tm, 2)


def calc_tm_difference(tm1, tm2):
    """
    Forward/Reverse primer의 Tm 차이(ΔTm)를 계산한다.
    - 일반적으로 ΔTm이 2℃ 이내일 때 두 프라이머가 비슷한 조건에서
      함께 결합(annealing)하기 유리하다고 참고된다.
    - 반환값: (ΔTm, 2℃ 이내 여부)
    """
    diff = round(abs(tm1 - tm2), 2)
    within_2c = diff <= 2.0
    return diff, within_2c


# ---------------------------------------------------------------------------
# 3. 결합 위치 탐색 (모든 위치) 및 PCR 산물 길이 계산
# ---------------------------------------------------------------------------

def find_all_binding_positions(target, primer, is_reverse=False):
    """
    표적 서열 내에서 프라이머가 결합 가능한 '모든' 위치를 찾는다.
    - is_reverse=False: forward primer를 그대로 표적 서열에서 탐색
    - is_reverse=True: reverse primer는 반대가닥에 결합하므로
      역상보서열(reverse complement)로 바꾼 뒤 탐색
    - str.find()를 반복 호출해서 겹치지 않는 모든 위치를 리스트로 반환
    - 하나도 찾지 못하면 빈 리스트를 반환한다.
    """
    if is_reverse:
        search_seq = str(Seq(primer).reverse_complement())
    else:
        search_seq = primer

    positions = []
    start = 0
    while True:
        pos = target.find(search_seq, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1  # 다음 탐색은 한 칸 뒤에서부터 (겹치는 위치도 찾기 위함)

    return positions, search_seq


def calc_pcr_product_length(target, forward_primer, reverse_primer):
    """
    PCR 산물 길이를 계산한다.

    개선 사항:
    - forward/reverse 각각 '모든' 결합 위치 후보를 구한 뒤,
      모든 조합(forward 위치 x reverse 위치)을 검토한다.
    - 올바른 PCR 증폭이 되려면 forward primer가 reverse primer보다
      앞쪽(더 작은 인덱스)에서 시작해야 하고, reverse primer의 결합
      부위(역상보서열) 끝이 forward 시작 위치보다 뒤에 있어야 한다.
      (= 두 프라이머가 서로 마주보며 가운데 영역을 증폭하는 배치)
    - 조건을 만족하는 조합이 하나도 없으면 산물 길이를 계산하지 않고
      그 이유를 함께 반환한다.
    - 조건을 만족하는 조합이 여러 개면, 가장 짧은 산물(가장 안쪽 조합)을
      대표 결과로 반환하고 나머지 조합도 함께 알려준다.
    """
    fwd_positions, fwd_seq = find_all_binding_positions(target, forward_primer, is_reverse=False)
    rev_positions, rev_seq = find_all_binding_positions(target, reverse_primer, is_reverse=True)

    result = {
        "fwd_positions": fwd_positions,
        "rev_positions": rev_positions,
        "valid_combinations": [],  # (fwd_pos, rev_pos, product_length)
        "product_length": None,
        "reason": "",
    }

    if not fwd_positions and not rev_positions:
        result["reason"] = "Forward와 Reverse primer 모두 표적 서열에서 결합 위치를 찾지 못했습니다."
        return result

    if not fwd_positions:
        result["reason"] = "Forward primer의 결합 위치를 찾지 못했습니다."
        return result

    if not rev_positions:
        result["reason"] = "Reverse primer의 결합 위치를 찾지 못했습니다."
        return result

    # 모든 (forward 위치, reverse 위치) 조합을 검토한다.
    for fwd_pos in fwd_positions:
        for rev_pos in rev_positions:
            rev_end = rev_pos + len(rev_seq)
            # 조건: forward가 reverse보다 앞에서 시작하고, 두 프라이머가 마주보는 배치여야 함
            if fwd_pos < rev_end and rev_pos >= fwd_pos:
                product_length = rev_end - fwd_pos
                result["valid_combinations"].append((fwd_pos, rev_pos, product_length))

    if not result["valid_combinations"]:
        result["reason"] = (
            "Forward와 Reverse primer의 결합 위치가 올바른 증폭 방향(서로 마주보는 배치)을 "
            "만족하지 않습니다. (Forward가 Reverse보다 뒤에 있거나 겹치는 등)"
        )
        return result

    # 가장 짧은(가장 안쪽) 조합을 대표 산물로 선택
    best = min(result["valid_combinations"], key=lambda x: x[2])
    result["product_length"] = best[2]
    result["best_combination"] = best
    return result


# ---------------------------------------------------------------------------
# 4. 구조적 특성 분석 (Hairpin, Self-dimer, Heterodimer)
# ---------------------------------------------------------------------------

def analyze_structure(primer, other_primer=None):
    """
    Primer3-py를 이용해 프라이머의 구조적 특성을 분석한다.
    - Hairpin: 프라이머 자체가 접혀서 생기는 2차 구조
    - Self-dimer(Homodimer): 같은 프라이머 두 분자끼리 결합하는 구조
    - Heterodimer: forward와 reverse primer가 서로 결합하는 구조
    각 구조는 ΔG(결합 자유에너지)와 Tm(해당 구조가 안정적으로 유지되는 온도)로 평가한다.
    ΔG가 음수로 클수록 구조가 더 안정적으로 형성된다는 뜻이며,
    이런 구조가 강하게 형성되면 프라이머가 표적 서열 대신 자기 자신이나
    서로 결합하는 데 쓰여 PCR 효율이 떨어질 가능성이 있다고 참고된다.
    (단, -9 kcal/mol은 절대적인 합격/불합격 기준이 아니라 참고용 경고 기준이다.)
    """
    hairpin = primer3.bindings.calc_hairpin(primer)
    homodimer = primer3.bindings.calc_homodimer(primer)

    result = {
        "hairpin_dg": round(hairpin.dg / 1000, 2),       # cal/mol -> kcal/mol
        "hairpin_tm": round(hairpin.tm, 2),
        "homodimer_dg": round(homodimer.dg / 1000, 2),
        "homodimer_tm": round(homodimer.tm, 2),
        "heterodimer_dg": None,
        "heterodimer_tm": None,
    }

    if other_primer is not None:
        heterodimer = primer3.bindings.calc_heterodimer(primer, other_primer)
        result["heterodimer_dg"] = round(heterodimer.dg / 1000, 2)
        result["heterodimer_tm"] = round(heterodimer.tm, 2)

    return result


# ---------------------------------------------------------------------------
# 5. 최종 결과 요약 및 기준 충족 여부 출력
# ---------------------------------------------------------------------------

DG_WARNING_THRESHOLD = -9.0  # 참고용 경고 기준 (절대 기준 아님)


def status_text(is_ok):
    return "양호" if is_ok else "주의"


def evaluate_primer(name, length, gc, tm, structure):
    """
    한 프라이머(forward 또는 reverse)의 길이/GC/Tm/구조 기준 충족 여부를 출력한다.
    (실험 성공을 보장하는 판정이 아니라, 설계 시 참고용 기준 충족 여부임)
    """
    print(f"\n--- {name} Primer ---")

    length_ok = 18 <= length <= 24
    print(f"길이: {length} bp  [{status_text(length_ok)}] (권장: 18~24bp)")

    gc_ok = 40 <= gc <= 60
    print(f"GC 함량: {gc} %  [{status_text(gc_ok)}] (권장: 40~60%)")

    tm_ok = 50 <= tm <= 65
    print(f"Tm (열역학적 계산): {tm} ℃  [{status_text(tm_ok)}] (권장: 50~65℃)")

    for key_dg, key_tm, label in [
        ("hairpin_dg", "hairpin_tm", "Hairpin"),
        ("homodimer_dg", "homodimer_tm", "Self-dimer"),
        ("heterodimer_dg", "heterodimer_tm", "Heterodimer"),
    ]:
        dg = structure[key_dg]
        tm_struct = structure[key_tm]
        if dg is None:
            continue
        dg_ok = dg > DG_WARNING_THRESHOLD
        print(
            f"{label}: ΔG {dg} kcal/mol, Tm {tm_struct} ℃  "
            f"[{status_text(dg_ok)}] (참고 경고 기준: {DG_WARNING_THRESHOLD} kcal/mol 이하면 "
            f"구조가 안정적으로 형성되어 PCR 효율에 영향을 줄 수 있음)"
        )


def print_final_summary(
    forward_primer, reverse_primer,
    fwd_len, rev_len, fwd_gc, rev_gc, fwd_tm, rev_tm,
    fwd_structure, rev_structure,
    product_result,
):
    """
    프라이머 분석 전체 결과를 한 번에 정리해서 출력한다.
    사용자가 어떤 부분이 양호하고 어떤 부분을 수정해야 하는지 한눈에 볼 수 있도록
    항목별로 상태([양호]/[주의])를 함께 표시한다.
    """
    print("\n================ 최종 결과 요약 ================")

    evaluate_primer("Forward", fwd_len, fwd_gc, fwd_tm, fwd_structure)
    evaluate_primer("Reverse", rev_len, rev_gc, rev_tm, rev_structure)

    # Tm 차이
    tm_diff, tm_diff_ok = calc_tm_difference(fwd_tm, rev_tm)
    print(f"\n--- Forward/Reverse Tm 차이 ---")
    print(f"ΔTm: {tm_diff} ℃  [{status_text(tm_diff_ok)}] (권장: 2℃ 이내)")
    if not tm_diff_ok:
        print("주의: Tm 차이가 커서 두 프라이머가 같은 annealing 온도 조건에서 "
              "동시에 잘 결합하지 않을 가능성이 있습니다.")

    # 결합 위치 및 산물 길이
    print(f"\n--- 결합 위치 및 PCR 산물 길이 ---")
    fwd_positions = product_result["fwd_positions"]
    rev_positions = product_result["rev_positions"]
    print(f"Forward primer 결합 위치: {fwd_positions if fwd_positions else '없음'}")
    print(f"Reverse primer 결합 위치(역상보서열 기준): {rev_positions if rev_positions else '없음'}")

    if product_result["product_length"] is not None:
        best_fwd, best_rev, best_len = product_result["best_combination"]
        print(f"[양호] 예상 PCR 산물 길이: {best_len} bp "
              f"(Forward 위치 {best_fwd}, Reverse 위치 {best_rev} 조합 기준)")
        if len(product_result["valid_combinations"]) > 1:
            print(f"참고: 조건을 만족하는 결합 위치 조합이 "
                  f"{len(product_result['valid_combinations'])}개 발견되었습니다. "
                  f"가장 짧은 산물을 기준으로 표시했습니다.")
            print(f"전체 조합: {product_result['valid_combinations']}")
    else:
        print(f"[주의] PCR 산물 길이를 계산할 수 없습니다. 이유: {product_result['reason']}")

    print("\n※ 위 기준들은 일반적으로 참고되는 설계 가이드라인이며, 실제 PCR 성공 여부를 보장하지 않습니다.")
    print("※ [주의]로 표시된 항목은 프라이머를 재설계할 때 우선적으로 검토해 볼 부분입니다.")


# ---------------------------------------------------------------------------
# 6. 메인 실행부
# ---------------------------------------------------------------------------

def main():
    target_dna = get_sequence_input("표적 DNA 서열을 입력하세요: ")
    forward_primer = get_sequence_input("Forward primer 서열을 입력하세요: ")
    reverse_primer = get_sequence_input("Reverse primer 서열을 입력하세요: ")

    if not (target_dna and forward_primer and reverse_primer):
        print("\n서열 입력에 문제가 있어 다음 단계로 진행할 수 없습니다.")
        return

    print("\n모든 서열이 정상적으로 입력되었습니다.")

    # 기본 특성 계산
    fwd_len = calc_length(forward_primer)
    rev_len = calc_length(reverse_primer)
    fwd_gc = calc_gc_content(forward_primer)
    rev_gc = calc_gc_content(reverse_primer)
    fwd_tm_simple = calc_tm_simple(forward_primer)
    rev_tm_simple = calc_tm_simple(reverse_primer)
    fwd_tm_thermo = calc_tm_thermo(forward_primer)
    rev_tm_thermo = calc_tm_thermo(reverse_primer)

    print(f"\n[Forward Primer: {forward_primer}]")
    print(f"길이: {fwd_len} bp / GC 함량: {fwd_gc} %")
    print(f"Tm (단순 계산식): {fwd_tm_simple} ℃ / Tm (열역학적 계산): {fwd_tm_thermo} ℃")

    print(f"\n[Reverse Primer: {reverse_primer}]")
    print(f"길이: {rev_len} bp / GC 함량: {rev_gc} %")
    print(f"Tm (단순 계산식): {rev_tm_simple} ℃ / Tm (열역학적 계산): {rev_tm_thermo} ℃")

    # 결합 위치 및 PCR 산물 길이 (모든 위치 탐색 + 방향 검증 포함)
    product_result = calc_pcr_product_length(target_dna, forward_primer, reverse_primer)

    print("\n[결합 위치 확인]")
    print(f"Forward primer 결합 위치: {product_result['fwd_positions'] if product_result['fwd_positions'] else '찾을 수 없음'}")
    print(f"Reverse primer 결합 위치(역상보서열 기준): {product_result['rev_positions'] if product_result['rev_positions'] else '찾을 수 없음'}")

    if product_result["product_length"] is not None:
        best_fwd, best_rev, best_len = product_result["best_combination"]
        print(f"\n예상 PCR 산물 길이: {best_len} bp (Forward 위치 {best_fwd}, Reverse 위치 {best_rev})")
    else:
        print(f"\nPCR 산물 길이를 계산할 수 없습니다. 이유: {product_result['reason']}")

    # 구조 분석
    fwd_structure = analyze_structure(forward_primer, other_primer=reverse_primer)
    rev_structure = analyze_structure(reverse_primer, other_primer=forward_primer)

    print("\n[Forward Primer 구조 분석]")
    print(f"Hairpin ΔG: {fwd_structure['hairpin_dg']} kcal/mol (Tm: {fwd_structure['hairpin_tm']} ℃)")
    print(f"Self-dimer ΔG: {fwd_structure['homodimer_dg']} kcal/mol (Tm: {fwd_structure['homodimer_tm']} ℃)")
    print(f"Heterodimer ΔG: {fwd_structure['heterodimer_dg']} kcal/mol (Tm: {fwd_structure['heterodimer_tm']} ℃)")

    print("\n[Reverse Primer 구조 분석]")
    print(f"Hairpin ΔG: {rev_structure['hairpin_dg']} kcal/mol (Tm: {rev_structure['hairpin_tm']} ℃)")
    print(f"Self-dimer ΔG: {rev_structure['homodimer_dg']} kcal/mol (Tm: {rev_structure['homodimer_tm']} ℃)")
    print(f"Heterodimer ΔG: {rev_structure['heterodimer_dg']} kcal/mol (Tm: {rev_structure['heterodimer_tm']} ℃)")

    # 최종 요약 (모든 항목 통합)
    print_final_summary(
        forward_primer, reverse_primer,
        fwd_len, rev_len, fwd_gc, rev_gc, fwd_tm_thermo, rev_tm_thermo,
        fwd_structure, rev_structure,
        product_result,
    )


if __name__ == "__main__":
    main()
