"""
primer_analysis.py

목적: 사용자가 입력한 표적 DNA 서열과 Forward/Reverse primer를 분석하는 프로그램
- 서열 유효성 검사
- 프라이머 길이 / GC 함량 / Tm 계산
- 표적 DNA 내 결합 위치 및 예상 PCR 산물 길이 계산
- Hairpin / Self-dimer / Heterodimer 구조 분석
- 일반적인 설계 기준과 비교한 결과 요약 출력

사용 라이브러리: Biopython, primer3-py
주의: 본 프로그램의 판정은 일반적으로 통용되는 설계 가이드라인과의 비교일 뿐이며,
      실제 PCR 실험의 성공 여부를 보장하지 않는다.
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


# ---------------------------------------------------------------------------
# 3. 결합 위치 및 PCR 산물 길이 계산
# ---------------------------------------------------------------------------

def find_binding_position(target, primer, is_reverse=False):
    """
    표적 서열 내에서 프라이머의 결합 위치를 찾는다.
    reverse primer는 역상보서열로 변환한 뒤 탐색한다.
    못 찾으면 위치는 -1로 반환한다.
    """
    if is_reverse:
        search_seq = str(Seq(primer).reverse_complement())
    else:
        search_seq = primer

    position = target.find(search_seq)
    return position, search_seq


def calc_pcr_product_length(target, forward_primer, reverse_primer):
    """
    PCR 산물 길이를 계산한다.
    forward 결합 시작 위치 ~ reverse 결합부위(역상보서열) 끝 위치까지의 거리.
    둘 중 하나라도 결합 위치를 못 찾으면 None을 반환한다.
    """
    fwd_pos, fwd_seq = find_binding_position(target, forward_primer, is_reverse=False)
    rev_pos, rev_seq = find_binding_position(target, reverse_primer, is_reverse=True)

    if fwd_pos == -1 or rev_pos == -1:
        return None, fwd_pos, rev_pos

    product_end = rev_pos + len(rev_seq)
    product_length = product_end - fwd_pos
    return product_length, fwd_pos, rev_pos


# ---------------------------------------------------------------------------
# 4. 구조적 특성 분석 (Hairpin, Self-dimer, Heterodimer)
# ---------------------------------------------------------------------------

def analyze_structure(primer, other_primer=None):
    """
    Primer3-py를 이용해 프라이머의 구조적 특성을 분석한다.
    ΔG(결합 자유에너지, kcal/mol)가 음수로 클수록 구조가 안정적으로 형성된다는 뜻이며,
    PCR 효율에 악영향을 줄 가능성이 있다고 참고적으로 해석한다.
    """
    hairpin = primer3.bindings.calc_hairpin(primer)
    homodimer = primer3.bindings.calc_homodimer(primer)

    result = {
        "hairpin_dg": round(hairpin.dg / 1000, 2),
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

def evaluate_primer(name, length, gc, tm, structure):
    """
    계산된 값들을 일반적으로 통용되는 설계 기준과 비교해 충족 여부를 표시한다.
    (실험 성공을 보장하는 판정이 아니라, 설계 시 참고용 기준 충족 여부임)
    """
    print(f"\n=== {name} Primer 최종 요약 ===")

    length_note = "(기준 충족: 18~24bp)" if 18 <= length <= 24 else "(기준 18~24bp 권장)"
    print(f"길이: {length} bp  {length_note}")

    gc_note = "(기준 충족: 40~60%)" if 40 <= gc <= 60 else "(기준 40~60% 권장)"
    print(f"GC 함량: {gc} %  {gc_note}")

    tm_note = "(기준 충족: 50~65℃)" if 50 <= tm <= 65 else "(기준 50~65℃ 권장)"
    print(f"Tm(열역학적 계산): {tm} ℃  {tm_note}")

    for key, label in [
        ("hairpin_dg", "Hairpin"),
        ("homodimer_dg", "Self-dimer"),
        ("heterodimer_dg", "Heterodimer"),
    ]:
        value = structure[key]
        if value is None:
            continue
        note = "(주의: 구조 형성 가능성 있음)" if value <= -9 else "(양호)"
        print(f"{label} ΔG: {value} kcal/mol  {note}")


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

    # 결합 위치 및 PCR 산물 길이
    fwd_pos, fwd_search_seq = find_binding_position(target_dna, forward_primer, is_reverse=False)
    rev_pos, rev_search_seq = find_binding_position(target_dna, reverse_primer, is_reverse=True)

    print("\n[결합 위치 확인]")
    if fwd_pos == -1:
        print("Forward primer: 표적 서열에서 결합 위치를 찾을 수 없습니다.")
    else:
        print(f"Forward primer 결합 위치: {fwd_pos} (0-based index)")

    if rev_pos == -1:
        print(f"Reverse primer: 표적 서열에서 결합 위치를 찾을 수 없습니다. (역상보서열 기준: {rev_search_seq})")
    else:
        print(f"Reverse primer 결합 위치(역상보서열 기준): {rev_pos} (0-based index)")

    product_length, _, _ = calc_pcr_product_length(target_dna, forward_primer, reverse_primer)
    if product_length is not None:
        print(f"\n예상 PCR 산물 길이: {product_length} bp")
    else:
        print("\nPCR 산물 길이를 계산할 수 없습니다. (프라이머 결합 위치 확인 실패)")

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

    # 최종 요약
    evaluate_primer("Forward", fwd_len, fwd_gc, fwd_tm_thermo, fwd_structure)
    evaluate_primer("Reverse", rev_len, rev_gc, rev_tm_thermo, rev_structure)

    print("\n※ 위 기준은 일반적으로 참고되는 설계 가이드라인이며, 실제 PCR 성공 여부를 보장하지 않습니다.")


if __name__ == "__main__":
    main()
