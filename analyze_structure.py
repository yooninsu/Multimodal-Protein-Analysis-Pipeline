#!/usr/bin/env python3
"""
AlphaFold 예측 구조와 에피토프 정보를 분석하는 스크립트
CSV 열 이름: uniprot_pos, uniprot_aa, label_epitope 지원
"""

import argparse
import pandas as pd
import numpy as np
from Bio import PDB
from Bio.PDB import SASA
import warnings
warnings.filterwarnings('ignore')

def calculate_sasa(structure):
    """
    각 잔기의 Solvent Accessible Surface Area 계산
    """
    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R")
    
    sasa_values = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # 표준 잔기만
                    res_num = residue.id[1]
                    sasa_values[res_num] = residue.sasa
    
    return sasa_values

def extract_plddt_values(structure):
    """
    B-factor 필드에서 pLDDT 값 추출
    """
    plddt_values = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':  # 표준 잔기만
                    res_num = residue.id[1]
                    # 각 잔기의 평균 B-factor (pLDDT) 계산
                    b_factors = [atom.bfactor for atom in residue]
                    plddt_values[res_num] = np.mean(b_factors)
    
    return plddt_values

def calculate_secondary_structure(structure):
    """
    간단한 2차 구조 분류
    """
    ss_values = {}
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == ' ':
                    res_num = residue.id[1]
                    ss_values[res_num] = 'C'  # Coil (기본값)
    
    return ss_values

def analyze_epitopes(pdb_file, truth_csv, output_file):
    """
    메인 분석 함수
    """
    print(f"Loading PDB file: {pdb_file}")
    
    # PDB 파일 파싱
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    # 구조적 특징 계산
    print("Calculating SASA values...")
    sasa_values = calculate_sasa(structure)
    
    print("Extracting pLDDT values...")
    plddt_values = extract_plddt_values(structure)
    
    print("Analyzing secondary structure...")
    ss_values = calculate_secondary_structure(structure)
    
    # 정답 CSV 읽기
    print(f"Loading truth CSV: {truth_csv}")
    truth_df = pd.read_csv(truth_csv)
    
    print("Available columns:", truth_df.columns.tolist())
    
    # 열 이름 매핑 (유연하게 처리)
    column_mapping = {}
    if 'uniprot_pos' in truth_df.columns:
        column_mapping['residue_index'] = 'uniprot_pos'
    elif 'residue_index' in truth_df.columns:
        column_mapping['residue_index'] = 'residue_index'
    
    if 'uniprot_aa' in truth_df.columns:
        column_mapping['residue_name'] = 'uniprot_aa'
    elif 'residue_name' in truth_df.columns:
        column_mapping['residue_name'] = 'residue_name'
    
    if 'label_epitope' in truth_df.columns:
        column_mapping['is_epitope'] = 'label_epitope'
    elif 'is_epitope' in truth_df.columns:
        column_mapping['is_epitope'] = 'is_epitope'
    
    print("Column mapping:", column_mapping)
    
    # 결과 데이터프레임 생성
    results = []
    
    for _, row in truth_df.iterrows():
        # 매핑된 열 이름 사용
        res_idx = row[column_mapping['residue_index']]
        
        # NaN 값 처리
        if pd.isna(res_idx):
            continue
            
        result_row = {
            'residue_index': int(res_idx) if not pd.isna(res_idx) else 0,
            'residue_name': row[column_mapping['residue_name']] if column_mapping['residue_name'] in row else 'UNK',
            'is_epitope': int(row[column_mapping['is_epitope']]) if column_mapping['is_epitope'] in row and not pd.isna(row[column_mapping['is_epitope']]) else 0,
            'sasa': sasa_values.get(int(res_idx) if not pd.isna(res_idx) else 0, 0),
            'plddt': plddt_values.get(int(res_idx) if not pd.isna(res_idx) else 0, 0),
            'secondary_structure': ss_values.get(int(res_idx) if not pd.isna(res_idx) else 0, 'Unknown'),
            'original_sasa_complex': row.get('sasa_complex', 0),
            'original_sasa_monomer': row.get('sasa_monomer', 0),
        }
        results.append(result_row)
    
    # 결과 정렬 및 저장
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('residue_index')
    
    # 통계 계산
    epitope_count = len(results_df[results_df['is_epitope'] == 1])
    total_count = len(results_df)
    
    print("\n=== Analysis Summary ===")
    print(f"Total residues analyzed: {total_count}")
    print(f"Epitope residues: {epitope_count}")
    print(f"Non-epitope residues: {total_count - epitope_count}")
    
    if epitope_count > 0:
        epitope_stats = results_df[results_df['is_epitope'] == 1]
        non_epitope_stats = results_df[results_df['is_epitope'] == 0]
        
        print(f"\nMean pLDDT (epitope): {epitope_stats['plddt'].mean():.2f}")
        print(f"Mean pLDDT (non-epitope): {non_epitope_stats['plddt'].mean():.2f}")
        
        print(f"\nMean SASA (epitope): {epitope_stats['sasa'].mean():.2f}")
        print(f"Mean SASA (non-epitope): {non_epitope_stats['sasa'].mean():.2f}")
    
    # 결과 저장
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # 요약 파일 생성
    summary_file = output_file.replace('.csv', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("=== AlphaFold Epitope Analysis Summary ===\n\n")
        f.write(f"Input PDB: {pdb_file}\n")
        f.write(f"Truth CSV: {truth_csv}\n\n")
        f.write(f"Total residues: {total_count}\n")
        f.write(f"Epitope residues: {epitope_count}\n")
        f.write(f"Non-epitope residues: {total_count - epitope_count}\n\n")
        
        if epitope_count > 0:
            f.write(f"Mean pLDDT (epitope): {epitope_stats['plddt'].mean():.2f}\n")
            f.write(f"Mean pLDDT (non-epitope): {non_epitope_stats['plddt'].mean():.2f}\n")
            f.write(f"Mean SASA (epitope): {epitope_stats['sasa'].mean():.2f}\n")
            f.write(f"Mean SASA (non-epitope): {non_epitope_stats['sasa'].mean():.2f}\n")
    
    print(f"Summary saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Analyze AlphaFold structure for epitopes')
    parser.add_argument('--pdb_file', required=True, help='Path to AlphaFold PDB file')
    parser.add_argument('--truth_file', required=True, help='Path to epitope truth CSV')
    parser.add_argument('--output_file', required=True, help='Output CSV file path')
    
    args = parser.parse_args()
    
    analyze_epitopes(args.pdb_file, args.truth_file, args.output_file)

if __name__ == "__main__":
    main()