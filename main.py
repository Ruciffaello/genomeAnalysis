import tkinter as tk
import threading
from tkinter import filedialog
from tkinter import scrolledtext
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

from Bio import SeqIO

class RecordInfo:
    def __init__(self, record):
        self.id = record.id
        self.seq = record.seq
        self.length = len(record.seq)

class App:
    def __init__(self,window):
        self.window = window
        self.window.title('Genome Sequence Compare')
        self.records =[]

        #variable
        self.match_score = tk.StringVar()
        self.match_score.set("1.0")
        self.mismatch_score = tk.StringVar()
        self.mismatch_score.set("-1.0")
        self.open_gap_score = tk.StringVar()
        self.open_gap_score.set("-5.0")
        self.extend_gap_score = tk.StringVar()
        self.extend_gap_score.set("-0.5")

        # 建立兩個標籤，分別顯示不同按鈕選的檔案

        # ===== 主框架設定 =====
        top_frame = tk.Frame(window)
        top_frame.grid(row=0, column=0, padx=10, pady=10, columnspan=2, sticky="nsew")

        middle_frame = tk.Frame(window)
        middle_frame.grid(row=1, column=0, padx=10, pady=5, columnspan=2, sticky="nsew")

        bottom_frame = tk.Frame(window)
        bottom_frame.grid(row=2, column=0, padx=10, pady=10, columnspan=2, sticky="ew")

        # ===== 上半部 - 左側垂直兩個 Text（合併 text_output1+2） =====
        left_text_frame = tk.Frame(top_frame)  # ⭐ 新增這個 frame 包兩個 text
        left_text_frame.grid(row=0, column=0, rowspan=1, sticky="nsew")  # 放在原本 text_output1 的位置
        top_frame.rowconfigure(0, weight=1)
        top_frame.columnconfigure(0, weight=1)

        # 第一個 Text
        text_output1 = tk.Text(left_text_frame, height=11, width=90)  # 注意高度減半
        scroll1 = tk.Scrollbar(left_text_frame, command=text_output1.yview)
        text_output1.configure(yscrollcommand=scroll1.set)
        text_output1.grid(row=0, column=0, padx=5, sticky="nsew")
        scroll1.grid(row=0, column=1, sticky="ns")

        # 第二個 Text
        text_output2 = tk.Text(left_text_frame, height=11, width=90)  # 注意高度減半
        scroll2 = tk.Scrollbar(left_text_frame, command=text_output2.yview)
        text_output2.configure(yscrollcommand=scroll2.set)
        text_output2.grid(row=1, column=0, padx=5, sticky="nsew")
        scroll2.grid(row=1, column=1, sticky="ns")

        # 讓上下平均拉伸
        left_text_frame.rowconfigure(0, weight=1)
        left_text_frame.rowconfigure(1, weight=1)
        left_text_frame.columnconfigure(0, weight=1)

        # ===== 右側單一 Text3（維持不動） =====
        text_output3 = tk.Text(top_frame, height=22, width=100)
        scroll3 = tk.Scrollbar(top_frame, command=text_output3.yview)
        text_output3.configure(yscrollcommand=scroll3.set)
        text_output3.grid(row=0, column=2, padx=5, sticky="nsew")
        scroll3.grid(row=0, column=3, sticky="ns")

        # ===== 底部 - 檔案按鈕、參數輸入區、比對按鈕 =====
        # 檔案區塊
        file_frame1 = tk.Frame(bottom_frame)
        file_frame1.grid(row=0, column=0, padx=10, pady=5)

        file_label1 = tk.Label(file_frame1, text="檔案1尚未選取")
        file_label1.pack(side=tk.LEFT)
        file_btn1 = tk.Button(file_frame1, text="選擇檔案1", command=lambda: self.select_file(file_label1, text_output1))
        file_btn1.pack(side=tk.LEFT)

        file_frame2 = tk.Frame(bottom_frame)
        file_frame2.grid(row=1, column=0, padx=10, pady=5)

        file_label2 = tk.Label(file_frame2, text="檔案2尚未選取")
        file_label2.pack(side=tk.LEFT)
        file_btn2 = tk.Button(file_frame2, text="選擇檔案2", command=lambda: self.select_file(file_label2, text_output2))
        file_btn2.pack(side=tk.LEFT)

        # 對齊參數欄位
        param_frame = tk.LabelFrame(bottom_frame, text="比對參數", padx=10, pady=10)
        param_frame.grid(row=2, column=0, columnspan=2, sticky="w")

        labels = ["match_score", "mismatch_score", "open_gap_score", "extend_gap_score"]
        vars_ = [self.match_score, self.mismatch_score, self.open_gap_score, self.extend_gap_score]

        for i, (label, var) in enumerate(zip(labels, vars_)):
            tk.Label(param_frame, text=label, anchor="w", width=15).grid(row=i, column=0, sticky="w", pady=2)
            tk.Entry(param_frame, textvariable=var, width=20).grid(row=i, column=1, sticky="w", pady=2)

        # 比對按鈕
        self.compare_btn = tk.Button(bottom_frame, text="比對序列", command=lambda :self.start_alignment(self.records,text_output3,self.compare_btn))
        self.compare_btn.grid(row=4, column=1, sticky="e", pady=10)



    def select_file(self,target_label,target_widget):
        file_path = filedialog.askopenfilename(
            title="選擇檔案",
            filetypes=(("所有檔案", "*.*"), ("文字檔", "*.txt"), ("Fasta 檔案", "*.fasta"))

        )
        if file_path:
            file_name = file_path.split("/")[-1]  # 也可以用 os.path.basename(file_path)
            target_label.config(text=f"已選取: {file_name}")
            target_widget.delete("1.0", tk.END)
            self.analysis_festa(file_path,target_widget)

    def analysis_festa(self,filepath,target_widget):
        for record in SeqIO.parse(filepath, "fasta"):
            target_widget.insert(tk.END, f"ID: {record.id}\n")
            target_widget.insert(tk.END, f"Sequence: {record.seq}\n")
            target_widget.insert(tk.END, f"Length: {len(record.seq)}\n")
            self.records.append(record.seq)

    def start_alignment(self, seq_pair, target_widget,compare_button, alignment_type="global"):
        compare_button.config(state=tk.DISABLED)
        target_widget.insert(tk.END, "分析中...\n")
        target_widget.see(tk.END)
        target_widget.update_idletasks()

        if len(seq_pair) < 2 :
            target_widget.insert(tk.END, "請確認基因序列是否都已經輸入...\n")
            return

        threading.Thread(target= lambda : self.threaded_alignment(seq_pair, target_widget, alignment_type)).start()

    def analyze_gene_similarity(self,gene_seq_str, target_widget, alignment_type="global"):
        gene1_seq = Seq(gene_seq_str[0].upper())
        gene2_seq = Seq(gene_seq_str[1].upper())

        aligner = PairwiseAligner()
        # aligner.match_score = 1.0
        # aligner.mismatch_score = -1.0
        # aligner.open_gap_score = -5.0
        # aligner.extend_gap_score = -0.5

        aligner.match_score = float(self.match_score.get())
        aligner.mismatch_score = float(self.mismatch_score.get())
        aligner.open_gap_score = -float(self.open_gap_score.get())
        aligner.extend_gap_score = float(self.extend_gap_score.get())

        result = ""

        if alignment_type == "global":
            aligner.mode = "global"
            result += f"執行全域比對 (Needleman-Wunsch)..."
        elif alignment_type == "local":
            aligner.mode = "local"
            result += f"執行局部比對 (Smith-Waterman)..."
        else:
            return "alignment_type 必須是 'global' 或 'local'"

        alignments = aligner.align(gene1_seq, gene2_seq)

        result += "\n--- 比對結果 ---"
        if not alignments:
            result += "沒有找到有效比對。"
            return result

        best_alignment = alignments[0]
        result += f"最佳比對分數: {best_alignment.score:.2f}"
        result += "\n比對後的序列：\n"
        result += str(best_alignment) + '\n'

        aligned_seq1 = str(best_alignment[0])
        aligned_seq2 = str(best_alignment[1])

        matches = 0
        aligned_length_without_gaps = 0
        for char1, char2 in zip(aligned_seq1, aligned_seq2):
            if char1 != '-' and char2 != '-':
                aligned_length_without_gaps += 1
                if char1 == char2:
                    matches += 1

        similarity_percentage = (matches / aligned_length_without_gaps) * 100 if aligned_length_without_gaps > 0 else 0

        result += f"匹配的核苷酸數量: {matches}\n"
        result += f"比對後有效長度: {aligned_length_without_gaps}\n"
        result += f"相似度百分比: {similarity_percentage:.2f}%\n"
        result += "\n--- 詳細比對信息 ---\n"
        result += f"序列1長度: {len(gene1_seq)}\n"
        result += f"序列2長度: {len(gene2_seq)}\n"

        return result


    def threaded_alignment(self, seq_pair, target_widget, alignment_type):
        result = self.analyze_gene_similarity(seq_pair, alignment_type)
        # 回到主執行緒更新 GUI
        self.window.after(0, lambda: self.display_alignment_result(result, target_widget))

    def display_alignment_result(self, result_str, target_widget):
        target_widget.insert(tk.END, result_str + "\n")
        target_widget.see(tk.END)
        self.compare_btn.config(state=tk.NORMAL)



if __name__ == "__main__":
    window = tk.Tk()
    window.title('Biopython Toos')
    window.geometry('1400x550')

    window.columnconfigure(0, weight=1)  # 左邊區塊
    window.columnconfigure(1, weight=1)  # 右邊區塊
    window.rowconfigure(0, weight=1)

    # 啟動 GUI 事件迴圈
    app = App(window)
    window.mainloop()





