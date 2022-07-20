import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
from Prop3D.util import safe_remove
from Prop3D.parsers.container import Container

def parallel_align_seq2aln(obj, sequence, alignment):
    return obj.align_sequence_to_alignment(sequence, alignment, as_biopython=True)

class Muscle(Container):
    IMAGE = 'docker://edraizen/muscle:latest'
    LOCAL = ["muscle"]
    PARAMETERS = [
        (":in_file", "path:in", "in"),
        (":profile", "store_true"),
        (":maxiters", "str"),
        (":in1", "path:in"),
        (":in2", "path:in"),
        (":out_file", "path:out", "out"),
        (":quiet", "store_true"),
        (":log", "path:out"),
        (":diags", "store_true"),
        (":maxhours", "str"),
        (":version", "store_true")
    ]
    RETURN_FILES=True
    ARG_START="-"

    def align_sequence_to_alignment(self, sequence, alignment, as_biopython=False):
        if isinstance(sequence, SeqRecord):
            seq_file = os.path.join(self.work_dir, f"{sequence.id.replace('/', '_').replace('|', '-')}.fasta")
            with open(seq_file, "w") as f:
                SeqIO.write(sequence, f, "fasta")
            remove = True
        elif isinstance(seuence, str) and os.path.isfile(sequence):
            #Assume fasta
            seq_file = sequence
            remove = False
        else:
            raise RuntimeError("Invalid sequence")

        new_seq_aln = seq_file+".aln"
        self(profile=True, in1=alignment, in2=seq_file, out_file=new_seq_aln)
        for aligned_seq in SeqIO.parse(new_seq_aln, "fasta"):
            if aligned_seq.id == seq.id:
                break
        else:
            aligned_seq = None

        if aligned_seq is None:
            seq = next(SeqIO.parse(seq_file, "fasta"))
            aln = next(SeqIO.parse(alignment, "fasta"))
            aligned_seq = SeqRecord(Seq('-'*len(aln)), id=seq.id)

        if remove:
            safe_remove(seq_file)

        print("SUCCCES", aligned_seq.id)

        if as_biopython:
            return aligned_seq

        aligned_seq_file = os.path.join(self.work_dir, f"{aligned_seq.id.replace('/', '_')}.aln.fasta")
        with open(aligned_seq_file, "w") as f:
            SeqIO.write()

        return aligned_seq_file

    def align_sequences_to_alignment(self, sequences, alignment, n_jobs=1):
        new_alignment = os.path.join(self.work_dir, f"{os.path.splitext(os.path.basename(sequences))[0]}__{os.path.splitext(os.path.basename(alignment))[0]}.aln.fatsa")
        with open(new_alignment, "w") as new:
            for aligned_seq in Parallel(n_jobs=n_jobs)(delayed(parallel_align_seq2aln)(self, seq, alignment) for seq in SeqIO.parse(sequences, "fasta")):
                SeqIO.write(aligned_seq, new, "fasta")
        return new_alignment
