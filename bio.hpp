// C++ implementations of the functions bio specifies.
//
// A sourced header is included before the compiler's own declarations, so
// every type these functions speak is declared here and mapped in main.loc,
// and no runtime type (such as mlc::Unit) can be named.
//
// Other modules include this file as "bio-cpp/bio.hpp" (it is exposed by
// package.yaml), while pools that source bio-cpp include its plane copy.
// Those are two files with the same contents, so the guard is a macro:
// #pragma once compares files, and would define every type twice.

#ifndef MORLOC_BIO_CPP_BIO_HPP
#define MORLOC_BIO_CPP_BIO_HPP

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <limits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <zlib.h>


// Types
// -----

struct BioGenome {};
struct BioTranscript {};

enum class BioStrand : uint8_t { Plus = 0, Minus = 1, Unstranded = 2, Unknown = 3 };
enum class BioOrient : uint8_t { Fwd = 0, Rev = 1 };

struct BioSegment {
    uint64_t start;
    uint64_t end;
    BioStrand strand;
};

template <class K> struct BioLoc {
    int seqid;
    BioSegment segment;
};

struct BioBlock {
    uint64_t srcStart;
    uint64_t dstStart;
    uint64_t len;
    BioOrient orient;
};

template <class A, class B> struct BioMapping {
    int src;
    int dst;
    std::vector<BioBlock> blocks;
};

template <class K> struct BioPiece_Mapped;
template <class K> struct BioPiece_Unmapped;

template <class K> struct BioPiece {
    std::variant<std::shared_ptr<const BioPiece_Mapped<K>>,
                 std::shared_ptr<const BioPiece_Unmapped<K>>> v;
};

template <class K> struct BioPiece_Mapped { BioLoc<K> f0; };
template <class K> struct BioPiece_Unmapped { uint64_t f0; uint64_t f1; };

template <class H, class S> struct BioSequence {
    H header;
    S bases;
};

struct BioBaseQual {
    std::vector<uint8_t> bases;
    std::vector<uint8_t> quals;
};

enum class BioPhase : uint8_t { Phase0 = 0, Phase1 = 1, Phase2 = 2, PhaseNone = 3 };

using BioSpan = std::tuple<uint64_t, uint64_t>;

template <class N> struct BioFeature {
    N node;
    std::vector<BioSegment> parts;
    std::vector<BioFeature<N>> children;
};

template <class S, class N> struct BioRoot {
    S seqid;
    BioFeature<N> feature;
};

using BioAttrs = std::vector<std::tuple<std::string, std::vector<std::string>>>;

struct BioGffNode {
    std::string kind;
    std::string origin;
    std::string id;
    std::string attrs;
    std::vector<std::optional<double>> scores;
    std::vector<BioPhase> phases;
    bool created;
};

struct BioBed3 {
    std::string chrom;
    uint64_t start;
    uint64_t end;
};

struct BioBed6 {
    std::string chrom;
    uint64_t start;
    uint64_t end;
    std::string name;
    double score;
    BioStrand strand;
};

struct BioBed12 {
    std::string chrom;
    uint64_t start;
    uint64_t end;
    std::string name;
    double score;
    BioStrand strand;
    uint64_t thickStart;
    uint64_t thickEnd;
    std::tuple<bool, uint8_t, uint8_t, uint8_t> rgb;
    std::vector<std::tuple<uint64_t, uint64_t>> blocks;
};

using BioFastaEntry = BioSequence<std::string, std::vector<uint8_t>>;
using BioFastqEntry = BioSequence<std::string, BioBaseQual>;
using BioGffRoot = BioRoot<std::string, BioGffNode>;

struct BioCds {
    std::string id;
    BioSpan range;
    BioPhase phase;
};

struct BioGeneTx {
    std::string id;
    std::string kind;
    std::string origin;
    std::string attrs;
    std::vector<BioSpan> exons;
    std::optional<BioCds> cds;
};

struct BioGene {
    std::string seqid;
    std::string id;
    std::string kind;
    std::string origin;
    std::string attrs;
    BioSpan span;
    BioStrand strand;
    std::vector<BioGeneTx> transcripts;
};

struct BioComplexity {
    bool mixedStrand;
    bool splitFeature;
    bool missingExons;
    bool overlappingExons;
    bool spanMismatch;
    bool multipleCds;
    bool cdsNotInExons;
    bool phaseShift;
    bool underivedLines;
    bool extraChildren;
    bool extraColumns;
    bool anonymousTranscript;
};

struct BioGffEntry_SimpleCodingGene;
struct BioGffEntry_SimpleNonCodingGene;
struct BioGffEntry_ComplexCodingGene;
struct BioGffEntry_ComplexNonCodingGene;
struct BioGffEntry_OtherFeature;

struct BioGffEntry {
    std::variant<std::shared_ptr<const BioGffEntry_SimpleCodingGene>,
                 std::shared_ptr<const BioGffEntry_SimpleNonCodingGene>,
                 std::shared_ptr<const BioGffEntry_ComplexCodingGene>,
                 std::shared_ptr<const BioGffEntry_ComplexNonCodingGene>,
                 std::shared_ptr<const BioGffEntry_OtherFeature>> v;
};

struct BioGffEntry_SimpleCodingGene { BioGene f0; };
struct BioGffEntry_SimpleNonCodingGene { BioGene f0; };
struct BioGffEntry_ComplexCodingGene { BioComplexity f0; BioGffRoot f1; };
struct BioGffEntry_ComplexNonCodingGene { BioComplexity f0; BioGffRoot f1; };
struct BioGffEntry_OtherFeature { BioGffRoot f0; };


// Callback plumbing
// -----------------

// The result and first argument of a std::function type. A writer builds the
// sink it hands to a producer from the producer's own parameter type, which
// is how it learns the sink's result type without naming it.
template <class F> struct bio_fn;
template <class R, class... A> struct bio_fn<std::function<R(A...)>> {
    using ret = R;
    using arg0 = std::tuple_element_t<0, std::tuple<A...>>;
};

// Hand `f` to `produce` as its sink.
template <class Producer, class F>
void bio_drain(Producer& produce, F f) {
    using SinkFn = typename bio_fn<Producer>::arg0;
    using Ret = typename bio_fn<SinkFn>::ret;
    produce(SinkFn([&f](const auto& batch) {
        f(batch);
        return Ret{};
    }));
}

// Collects records and hands them to a sink in batches.
template <class T, class Sink>
class BioBatcher {
  public:
    BioBatcher(Sink& sink, size_t limit) : sink_(sink), limit_(limit) {}
    void add(T x, size_t weight) {
        weight_ += weight;
        batch_.push_back(std::move(x));
        if (weight_ >= limit_) flush();
    }
    void flush() {
        if (batch_.empty()) return;
        sink_(std::move(batch_));
        batch_.clear();
        weight_ = 0;
    }

  private:
    Sink& sink_;
    size_t limit_;
    size_t weight_ = 0;
    std::vector<T> batch_;
};

// Sequence batches hold about this many bytes; annotation batches this many
// trees.
constexpr size_t BIO_BATCH_BYTES = 1u << 20;
constexpr size_t BIO_BATCH_TREES = 4096;


// Coordinate mappings
// -------------------

inline BioStrand bio_flip(BioStrand s) {
    if (s == BioStrand::Plus) return BioStrand::Minus;
    if (s == BioStrand::Minus) return BioStrand::Plus;
    return s;
}

inline BioOrient bio_orient_mul(BioOrient a, BioOrient b) {
    return a == b ? BioOrient::Fwd : BioOrient::Rev;
}

inline void bio_sort_blocks(std::vector<BioBlock>& bs) {
    std::stable_sort(bs.begin(), bs.end(), [](const BioBlock& x, const BioBlock& y) {
        return x.srcStart < y.srcStart;
    });
}

// The image in B of A's [s, e), which must lie within the block.
inline std::pair<uint64_t, uint64_t> bio_block_image(const BioBlock& b, uint64_t s, uint64_t e) {
    if (b.orient == BioOrient::Fwd) {
        return {b.dstStart + (s - b.srcStart), b.dstStart + (e - b.srcStart)};
    }
    uint64_t srcEnd = b.srcStart + b.len;
    return {b.dstStart + (srcEnd - e), b.dstStart + (srcEnd - s)};
}

// A mapping whose kinds are fixed by the variable it initializes: the pool
// does not pass template arguments that appear only in a return type.
struct BioMappingAny {
    int src;
    int dst;
    std::vector<BioBlock> blocks;
    template <class A, class B> operator BioMapping<A, B>() const {
        return BioMapping<A, B>{src, dst, blocks};
    }
};

inline BioMappingAny bio_tiled(int src, int dst,
                               const std::vector<std::tuple<uint64_t, uint64_t, BioOrient>>& parts) {
    BioMappingAny m{src, dst, {}};
    uint64_t at = 0;
    for (const auto& [dstStart, len, orient] : parts) {
        m.blocks.push_back(BioBlock{at, dstStart, len, orient});
        at += len;
    }
    return m;
}

template <class A, class B>
BioMapping<B, A> bio_invert(const BioMapping<A, B>& m) {
    BioMapping<B, A> r{m.dst, m.src, {}};
    for (const auto& b : m.blocks) {
        r.blocks.push_back(BioBlock{b.dstStart, b.srcStart, b.len, b.orient});
    }
    bio_sort_blocks(r.blocks);
    return r;
}

template <class A, class B, class C>
BioMapping<A, C> bio_compose_blocks(const BioMapping<A, B>& ab, const BioMapping<B, C>& bc) {
    BioMapping<A, C> r{ab.src, bc.dst, {}};
    for (const auto& x : ab.blocks) {
        const uint64_t x0 = x.dstStart, x1 = x.dstStart + x.len;
        const BioBlock xInv{x.dstStart, x.srcStart, x.len, x.orient};
        for (const auto& y : bc.blocks) {
            const uint64_t is = std::max(x0, y.srcStart);
            const uint64_t ie = std::min(x1, y.srcStart + y.len);
            if (is >= ie) continue;
            // [is, ie) of B is the image under x of an A range and has an
            // image under y in C; the new block pairs the two starts.
            const auto a = bio_block_image(xInv, is, ie);
            const auto c = bio_block_image(y, is, ie);
            r.blocks.push_back(BioBlock{a.first, c.first, ie - is, bio_orient_mul(x.orient, y.orient)});
        }
    }
    bio_sort_blocks(r.blocks);
    return r;
}

// `project`, scanning blocks from index `first`; every block before it must
// end before the location (at or before a point location's position).
template <class A, class B>
std::vector<BioPiece<B>> bio_project_from(const BioMapping<A, B>& m, const BioLoc<A>& loc, size_t first) {
    if (loc.seqid != m.src) {
        throw std::runtime_error("project: the location is not on the mapping's source sequence");
    }
    const uint64_t s = loc.segment.start, e = loc.segment.end;
    if (s > e) throw std::runtime_error("project: segment start is past its end");

    std::vector<BioPiece<B>> out;
    auto mapped = [&](std::pair<uint64_t, uint64_t> img, BioOrient o) {
        const BioStrand st = o == BioOrient::Rev ? bio_flip(loc.segment.strand) : loc.segment.strand;
        const BioLoc<B> l{m.dst, BioSegment{img.first, img.second, st}};
        out.push_back(BioPiece<B>{std::make_shared<BioPiece_Mapped<B>>(BioPiece_Mapped<B>{l})});
    };
    auto unmapped = [&](uint64_t us, uint64_t ue) {
        out.push_back(BioPiece<B>{std::make_shared<BioPiece_Unmapped<B>>(BioPiece_Unmapped<B>{us, ue})});
    };

    if (s == e) {
        // A point between bases: the first block containing it, ends included.
        for (size_t i = first; i < m.blocks.size() && m.blocks[i].srcStart <= s; ++i) {
            const auto& b = m.blocks[i];
            if (s <= b.srcStart + b.len) {
                mapped(bio_block_image(b, s, s), b.orient);
                return out;
            }
        }
        unmapped(s, s);
        return out;
    }

    // Blocks are sorted by srcStart and may overlap. `covered` is the end of
    // the part of [s, e) some block has already covered.
    uint64_t covered = s;
    for (size_t i = first; i < m.blocks.size() && m.blocks[i].srcStart < e; ++i) {
        const auto& b = m.blocks[i];
        const uint64_t bs = std::max(s, b.srcStart);
        const uint64_t be = std::min(e, b.srcStart + b.len);
        if (bs >= be) continue;
        if (bs > covered) unmapped(covered, bs);
        mapped(bio_block_image(b, bs, be), b.orient);
        covered = std::max(covered, be);
    }
    if (covered < e) unmapped(covered, e);
    return out;
}

template <class A, class B>
std::vector<BioPiece<B>> bio_project(const BioMapping<A, B>& m, const BioLoc<A>& loc) {
    return bio_project_from(m, loc, 0);
}

// `project` of many locations through one mapping. The running maximum of
// block ends finds each location's first candidate block by binary search.
template <class A, class B>
std::vector<std::vector<BioPiece<B>>> bio_project_many(const BioMapping<A, B>& m, const std::vector<BioLoc<A>>& locs) {
    std::vector<uint64_t> reach(m.blocks.size());
    uint64_t hi = 0;
    for (size_t i = 0; i < m.blocks.size(); ++i) {
        if (i > 0 && m.blocks[i].srcStart < m.blocks[i - 1].srcStart)
            throw std::runtime_error("projectMany: the mapping's blocks are not sorted by srcStart");
        hi = std::max(hi, m.blocks[i].srcStart + m.blocks[i].len);
        reach[i] = hi;
    }
    std::vector<std::vector<BioPiece<B>>> out;
    out.reserve(locs.size());
    for (const auto& loc : locs) {
        const uint64_t s = loc.segment.start;
        const bool point = s == loc.segment.end;
        const auto it = std::partition_point(reach.begin(), reach.end(),
                                             [&](uint64_t r) { return point ? r < s : r <= s; });
        out.push_back(bio_project_from(m, loc, static_cast<size_t>(it - reach.begin())));
    }
    return out;
}

// The IUPAC complement of a base, case kept. Other bytes are unchanged.
inline uint8_t bio_complement(uint8_t c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A'; case 'U': return 'A';
        case 'G': return 'C'; case 'C': return 'G';
        case 'R': return 'Y'; case 'Y': return 'R';
        case 'K': return 'M'; case 'M': return 'K';
        case 'B': return 'V'; case 'V': return 'B';
        case 'D': return 'H'; case 'H': return 'D';
        case 'a': return 't'; case 't': return 'a'; case 'u': return 'a';
        case 'g': return 'c'; case 'c': return 'g';
        case 'r': return 'y'; case 'y': return 'r';
        case 'k': return 'm'; case 'm': return 'k';
        case 'b': return 'v'; case 'v': return 'b';
        case 'd': return 'h'; case 'h': return 'd';
        default: return c;
    }
}

template <class A, class B>
std::vector<uint8_t> bio_extract(const BioMapping<A, B>& m, const std::vector<uint8_t>& dst) {
    std::vector<uint8_t> out;
    for (const auto& b : m.blocks) {
        if (b.dstStart + b.len > dst.size()) {
            throw std::runtime_error("extract: a block runs past the end of the sequence");
        }
        const auto first = dst.begin() + b.dstStart;
        const auto last = first + b.len;
        if (b.orient == BioOrient::Fwd) {
            out.insert(out.end(), first, last);
        } else {
            for (auto it = last; it != first;) {
                --it;
                out.push_back(bio_complement(*it));
            }
        }
    }
    return out;
}


template <class A, class B>
std::vector<uint8_t> bio_extract_blocks(const BioMapping<A, B>& m, const std::vector<std::vector<uint8_t>>& slices) {
    if (slices.size() != m.blocks.size()) {
        throw std::runtime_error("extractBlocks: expected one slice per block");
    }
    std::vector<uint8_t> out;
    for (size_t i = 0; i < slices.size(); ++i) {
        const BioBlock& b = m.blocks[i];
        const auto& sl = slices[i];
        if (sl.size() != b.len) throw std::runtime_error("extractBlocks: a slice's length differs from its block's");
        if (b.orient == BioOrient::Fwd) {
            out.insert(out.end(), sl.begin(), sl.end());
        } else {
            for (auto it = sl.rbegin(); it != sl.rend(); ++it) out.push_back(bio_complement(*it));
        }
    }
    return out;
}


// Features
// --------

// Parts in transcript order: ascending, or descending when every part is on
// the minus strand. `perm` receives the new order of the old indices.
inline std::vector<size_t> bio_transcript_order(const std::vector<BioSegment>& parts) {
    std::vector<size_t> perm(parts.size());
    for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
    std::stable_sort(perm.begin(), perm.end(), [&](size_t a, size_t b) {
        return parts[a].start < parts[b].start;
    });
    const bool minus = !parts.empty() &&
        std::all_of(parts.begin(), parts.end(), [](const BioSegment& p) { return p.strand == BioStrand::Minus; });
    if (minus) std::reverse(perm.begin(), perm.end());
    return perm;
}

template <class T>
std::vector<T> bio_permute(const std::vector<T>& xs, const std::vector<size_t>& perm) {
    std::vector<T> out;
    out.reserve(xs.size());
    for (size_t i : perm) out.push_back(xs[i]);
    return out;
}

// A feature's parts as segments, in its part order.
template <class N>
const std::vector<BioSegment>& bio_segments(const BioFeature<N>& f) {
    return f.parts;
}

// A feature's parts as [start, end) spans, in its part order.
template <class N>
std::vector<BioSpan> bio_spans(const BioFeature<N>& f) {
    std::vector<BioSpan> out;
    out.reserve(f.parts.size());
    for (const auto& p : f.parts) out.emplace_back(p.start, p.end);
    return out;
}

// The strand every part shares, or nothing when they differ or there are none.
template <class N>
std::optional<BioStrand> bio_single_strand(const BioFeature<N>& f) {
    if (f.parts.empty()) return std::nullopt;
    for (const auto& p : f.parts)
        if (p.strand != f.parts[0].strand) return std::nullopt;
    return f.parts[0].strand;
}

// A feature from segments in any order, put in transcript order.
template <class N>
BioFeature<N> bio_feature(N node, const std::vector<BioSegment>& segs, std::vector<BioFeature<N>> children = {}) {
    return BioFeature<N>{std::move(node), bio_permute(segs, bio_transcript_order(segs)), std::move(children)};
}


// Views
// -----

struct BioTxModel {
    std::string id;
    std::string gene;
    std::string seqid;
    BioMapping<BioTranscript, BioGenome> mapping;
    std::optional<std::tuple<uint64_t, uint64_t>> cds;
};

// Collects the parts of `f`'s children of a kind.
template <class N>
std::vector<BioSegment> bio_child_parts(const BioFeature<N>& f, const std::string& kind) {
    std::vector<BioSegment> parts;
    for (const auto& c : f.children) {
        if (c.node.kind == kind) parts.insert(parts.end(), c.parts.begin(), c.parts.end());
    }
    return parts;
}

// A transcript model from its exons and CDS parts on genome sequence
// `dst`; the model's index in the output is `src`.
inline BioTxModel bio_tx_model(std::string id, std::string gene, std::string seqid, int src, int dst,
                               std::vector<BioSegment> exons, const std::vector<BioSegment>& cdsParts) {
    BioTxModel m;
    m.id = std::move(id);
    m.gene = std::move(gene);
    m.seqid = std::move(seqid);
    std::sort(exons.begin(), exons.end(), [](const BioSegment& a, const BioSegment& b) { return a.start < b.start; });
    const BioStrand st = exons[0].strand;
    for (const auto& e : exons) {
        if (e.strand != st || (st != BioStrand::Plus && st != BioStrand::Minus)) {
            throw std::runtime_error("transcriptModels: " + m.id + " has exons on mixed or unknown strands");
        }
    }
    if (st == BioStrand::Minus) std::reverse(exons.begin(), exons.end());
    const BioOrient o = st == BioStrand::Minus ? BioOrient::Rev : BioOrient::Fwd;
    std::vector<std::tuple<uint64_t, uint64_t, BioOrient>> parts;
    for (const auto& e : exons) parts.emplace_back(e.start, e.end - e.start, o);
    m.mapping = bio_tiled(src, dst, parts);
    if (!cdsParts.empty()) {
        const auto inv = bio_invert(m.mapping);
        uint64_t lo = UINT64_MAX, hi = 0;
        for (const auto& p : cdsParts) {
            const BioLoc<BioGenome> loc{dst, BioSegment{p.start, p.end, p.strand}};
            for (const auto& piece : bio_project(inv, loc)) {
                const auto* mp = std::get_if<std::shared_ptr<const BioPiece_Mapped<BioTranscript>>>(&piece.v);
                if (!mp) throw std::runtime_error("transcriptModels: a CDS of " + m.id + " lies outside its exons");
                lo = std::min(lo, (*mp)->f0.segment.start);
                hi = std::max(hi, (*mp)->f0.segment.end);
            }
        }
        m.cds = std::make_tuple(lo, hi);
    }
    return m;
}

inline int bio_seq_index(const std::unordered_map<std::string, int>& index, const std::string& seqid,
                         const std::string& what) {
    const auto it = index.find(seqid);
    if (it == index.end()) {
        throw std::runtime_error("transcriptModels: " + what + " is on " + seqid + ", which is not a genome sequence");
    }
    return it->second;
}

// The transcripts of one tree.
inline void bio_root_tx_models(const std::unordered_map<std::string, int>& index, const BioGffRoot& root,
                               std::vector<BioTxModel>& out) {
    const auto& g = root.feature;
    const std::string& gene = g.node.id;
    std::vector<const BioFeature<BioGffNode>*> txs;
    if (!bio_child_parts(g, "exon").empty()) txs.push_back(&g);
    for (const auto& c : g.children) {
        if (!bio_child_parts(c, "exon").empty()) txs.push_back(&c);
    }
    for (const BioFeature<BioGffNode>* t : txs) {
        const std::string id = !t->node.id.empty() ? t->node.id : gene + ":" + std::to_string(out.size());
        const int dst = bio_seq_index(index, root.seqid, id);
        out.push_back(bio_tx_model(id, gene, root.seqid, static_cast<int>(out.size()), dst,
                                   bio_child_parts(*t, "exon"), bio_child_parts(*t, "CDS")));
    }
}

// A standard transcript's CDS parts: its exons clipped to the coding range.
inline std::vector<BioSpan> bio_cds_spans(const BioGeneTx& t) {
    std::vector<BioSpan> out;
    if (!t.cds) return out;
    const auto [cs, ce] = t.cds->range;
    for (const auto& [s, e] : t.exons) {
        const uint64_t a = std::max(s, cs), b = std::min(e, ce);
        if (a < b) out.emplace_back(a, b);
    }
    return out;
}

template <class F>
void bio_entry_visit(const BioGffEntry& e, F onGene, std::function<void(const BioGffRoot&)> onRoot) {
    if (const auto* p = std::get_if<std::shared_ptr<const BioGffEntry_SimpleCodingGene>>(&e.v)) onGene((*p)->f0);
    else if (const auto* p = std::get_if<std::shared_ptr<const BioGffEntry_SimpleNonCodingGene>>(&e.v)) onGene((*p)->f0);
    else if (const auto* p = std::get_if<std::shared_ptr<const BioGffEntry_ComplexCodingGene>>(&e.v)) onRoot((*p)->f1);
    else if (const auto* p = std::get_if<std::shared_ptr<const BioGffEntry_ComplexNonCodingGene>>(&e.v)) onRoot((*p)->f1);
    else onRoot(std::get<std::shared_ptr<const BioGffEntry_OtherFeature>>(e.v)->f0);
}

inline std::string bio_entry_seqid(const BioGffEntry& e) {
    std::string s;
    bio_entry_visit(e, [&](const BioGene& g) { s = g.seqid; }, [&](const BioGffRoot& r) { s = r.seqid; });
    return s;
}

inline std::vector<BioTxModel> bio_transcript_models(const std::vector<std::string>& names,
                                                     const std::vector<BioGffEntry>& entries) {
    std::unordered_map<std::string, int> index;
    for (size_t i = 0; i < names.size(); ++i) index[names[i]] = static_cast<int>(i);
    std::vector<BioTxModel> out;
    for (const auto& e : entries) {
        bio_entry_visit(
            e,
            [&](const BioGene& g) {
                for (const auto& t : g.transcripts) {
                    const std::string id = !t.id.empty() ? t.id
                                           : !g.id.empty() ? g.id
                                                           : g.id + ":" + std::to_string(out.size());
                    const int dst = bio_seq_index(index, g.seqid, id);
                    std::vector<BioSegment> exons, cds;
                    for (const auto& [s, e2] : t.exons) exons.push_back(BioSegment{s, e2, g.strand});
                    for (const auto& [s, e2] : bio_cds_spans(t)) cds.push_back(BioSegment{s, e2, g.strand});
                    out.push_back(bio_tx_model(id, g.id, g.seqid, static_cast<int>(out.size()), dst, exons, cds));
                }
            },
            [&](const BioGffRoot& r) { bio_root_tx_models(index, r, out); });
    }
    return out;
}


// Text I/O
// --------

// Reads lines from a plain or gzip-compressed file (zlib reads both).
class BioLineReader {
  public:
    explicit BioLineReader(const std::string& path) : path_(path) {
        f_ = gzopen(path.c_str(), "rb");
        if (!f_) throw std::runtime_error("cannot open " + path);
        gzbuffer(f_, 1u << 18);
    }
    ~BioLineReader() {
        if (f_) gzclose(f_);
    }
    BioLineReader(const BioLineReader&) = delete;
    BioLineReader& operator=(const BioLineReader&) = delete;

    // The next line without its line ending; false at end of file.
    bool next(std::string& line) {
        line.clear();
        while (true) {
            if (pos_ == len_ && !fill()) {
                if (line.empty()) return false;
                if (line.back() == '\r') line.pop_back();
                ++lineno_;
                return true;
            }
            const char* start = buf_.data() + pos_;
            const void* nl = std::memchr(start, '\n', len_ - pos_);
            if (nl) {
                const size_t n = static_cast<const char*>(nl) - start;
                line.append(start, n);
                pos_ += n + 1;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                ++lineno_;
                return true;
            }
            line.append(start, len_ - pos_);
            pos_ = len_;
        }
    }

    // "path:line: msg", for the line last returned.
    std::runtime_error error(const std::string& msg) const {
        return std::runtime_error(path_ + ":" + std::to_string(lineno_) + ": " + msg);
    }

  private:
    bool fill() {
        const int n = gzread(f_, buf_.data(), static_cast<unsigned>(buf_.size()));
        if (n < 0) {
            int err;
            throw std::runtime_error("read error in " + path_ + ": " + gzerror(f_, &err));
        }
        pos_ = 0;
        len_ = static_cast<size_t>(n);
        return n > 0;
    }

    std::string path_;
    gzFile f_ = nullptr;
    std::vector<char> buf_ = std::vector<char>(1u << 18);
    size_t pos_ = 0, len_ = 0, lineno_ = 0;
};

// Writes a file, gzip-compressed when the path ends in ".gz". Output goes to
// "<path>.part" and is renamed to <path> by close(), so a write that fails
// part-way leaves no file that looks complete.
class BioWriter {
  public:
    explicit BioWriter(const std::string& path) : path_(path), part_(path + ".part") {
        gz_ = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz_) {
            g_ = gzopen(part_.c_str(), "wb");
        } else {
            p_ = std::fopen(part_.c_str(), "wb");
        }
        if (!g_ && !p_) throw std::runtime_error("cannot open " + part_ + " for writing");
    }
    ~BioWriter() {
        if (g_) gzclose(g_);
        if (p_) std::fclose(p_);
        if (!closed_) std::remove(part_.c_str());
    }
    BioWriter(const BioWriter&) = delete;
    BioWriter& operator=(const BioWriter&) = delete;

    void write(const char* data, size_t n) {
        if (n == 0) return;
        const bool ok = gz_ ? gzwrite(g_, data, static_cast<unsigned>(n)) == static_cast<int>(n)
                            : std::fwrite(data, 1, n, p_) == n;
        if (!ok) throw std::runtime_error("write error in " + path_);
    }
    void write(const std::string& s) { write(s.data(), s.size()); }
    void write(const std::vector<uint8_t>& v) {
        write(reinterpret_cast<const char*>(v.data()), v.size());
    }
    void put(char c) { write(&c, 1); }

    // Flush and close, reporting the errors a destructor would swallow.
    void close() {
        if (g_) {
            const int rc = gzclose(g_);
            g_ = nullptr;
            if (rc != Z_OK) throw std::runtime_error("error closing " + path_);
        }
        if (p_) {
            const int rc = std::fclose(p_);
            p_ = nullptr;
            if (rc != 0) throw std::runtime_error("error closing " + path_);
        }
        if (std::rename(part_.c_str(), path_.c_str()) != 0) {
            throw std::runtime_error("cannot rename " + part_ + " to " + path_);
        }
        closed_ = true;
    }

  private:
    std::string path_;
    std::string part_;
    bool closed_ = false;
    bool gz_ = false;
    gzFile g_ = nullptr;
    std::FILE* p_ = nullptr;
};

// Appends text to a file, as a new gzip member when the path ends in ".gz"
// (a reader of gzip reads concatenated members as one stream).
inline void bio_append_text(const std::string& path, const std::string& text) {
    const bool gz = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
    if (gz) {
        gzFile g = gzopen(path.c_str(), "ab");
        if (!g) throw std::runtime_error("cannot open " + path + " for appending");
        const bool ok = text.empty() || gzwrite(g, text.data(), static_cast<unsigned>(text.size())) == static_cast<int>(text.size());
        if (gzclose(g) != Z_OK || !ok) throw std::runtime_error("write error in " + path);
        return;
    }
    std::FILE* p = std::fopen(path.c_str(), "ab");
    if (!p) throw std::runtime_error("cannot open " + path + " for appending");
    const bool ok = std::fwrite(text.data(), 1, text.size(), p) == text.size();
    if (std::fclose(p) != 0 || !ok) throw std::runtime_error("write error in " + path);
}

// Creates an empty file, replacing any file at the path.
inline void bio_reset_file(const std::string& path) {
    std::FILE* p = std::fopen(path.c_str(), "wb");
    if (!p || std::fclose(p) != 0) throw std::runtime_error("cannot create " + path);
}

inline std::vector<std::string> bio_split(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t at = 0;
    while (true) {
        const size_t k = s.find(sep, at);
        if (k == std::string::npos) {
            out.push_back(s.substr(at));
            return out;
        }
        out.push_back(s.substr(at, k - at));
        at = k + 1;
    }
}

inline std::string bio_trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

// A non-negative integer field, or an error naming it.
inline uint64_t bio_parse_u64(const BioLineReader& in, const std::string& s, const char* what) {
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos) {
        throw in.error(std::string("bad ") + what + " '" + s + "'");
    }
    errno = 0;
    const unsigned long long v = std::strtoull(s.c_str(), nullptr, 10);
    if (errno == ERANGE) throw in.error(std::string(what) + " out of range '" + s + "'");
    return static_cast<uint64_t>(v);
}

// A number, or nothing for ".".
inline std::optional<double> bio_parse_score(const BioLineReader& in, const std::string& s) {
    if (s == ".") return std::nullopt;
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(s.c_str(), &end);
    if (s.empty() || *end != '\0' || errno == ERANGE) throw in.error("bad score '" + s + "'");
    return v;
}

// An integral value as an integer, else the shortest decimal text that
// reads back as the same double.
inline std::string bio_format_score(const std::optional<double>& v) {
    if (!v) return ".";
    char buf[64];
    if (*v == static_cast<double>(static_cast<long long>(*v)) && std::fabs(*v) < 1e15) {
        std::snprintf(buf, sizeof buf, "%lld", static_cast<long long>(*v));
        return buf;
    }
    for (int prec = 1; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof buf, "%.*g", prec, *v);
        if (std::strtod(buf, nullptr) == *v) break;
    }
    return buf;
}

inline BioStrand bio_parse_strand(const BioLineReader& in, const std::string& s) {
    if (s == "+") return BioStrand::Plus;
    if (s == "-") return BioStrand::Minus;
    if (s == ".") return BioStrand::Unstranded;
    if (s == "?") return BioStrand::Unknown;
    throw in.error("bad strand '" + s + "'");
}

inline char bio_strand_char(BioStrand s) {
    switch (s) {
        case BioStrand::Plus: return '+';
        case BioStrand::Minus: return '-';
        case BioStrand::Unstranded: return '.';
        default: return '?';
    }
}

// FASTA
// -----

template <class Sink>
void bio_read_fasta(const std::string& path, Sink sink) {
    BioLineReader in(path);
    BioBatcher<BioFastaEntry, Sink> out(sink, BIO_BATCH_BYTES);
    std::string line;
    bool open = false;
    BioFastaEntry cur;
    auto finish = [&]() {
        const size_t w = cur.header.size() + cur.bases.size();
        out.add(std::move(cur), w);
        cur = BioFastaEntry{};
    };
    while (in.next(line)) {
        if (line.empty() || line[0] == ';') continue;
        if (line[0] == '>') {
            if (open) finish();
            cur.header = line.substr(1);
            open = true;
            continue;
        }
        if (!open) throw in.error("sequence before the first '>' header");
        // memchr is vectorized in the C library
        if (!std::memchr(line.data(), ' ', line.size()) && !std::memchr(line.data(), '\t', line.size())) {
            cur.bases.insert(cur.bases.end(), line.begin(), line.end());
            continue;
        }
        for (char c : line) {
            if (c != ' ' && c != '\t') cur.bases.push_back(static_cast<uint8_t>(c));
        }
    }
    if (open) finish();
    out.flush();
}

// Sequence lines are wrapped at this width.
constexpr size_t BIO_FASTA_WIDTH = 60;

inline void bio_fasta_text(std::string& s, const BioFastaEntry& e) {
    s += '>';
    s += e.header;
    s += '\n';
    for (size_t i = 0; i < e.bases.size(); i += BIO_FASTA_WIDTH) {
        const size_t n = std::min(BIO_FASTA_WIDTH, e.bases.size() - i);
        s.append(reinterpret_cast<const char*>(e.bases.data()) + i, n);
        s += '\n';
    }
}

inline std::string bio_format_fasta(const std::vector<BioFastaEntry>& batch) {
    std::string s;
    for (const auto& e : batch) bio_fasta_text(s, e);
    return s;
}

template <class Producer>
void bio_write_fasta(const std::string& path, Producer produce) {
    BioWriter out(path);
    bio_drain(produce, [&](const std::vector<BioFastaEntry>& batch) { out.write(bio_format_fasta(batch)); });
    out.close();
}


// FASTQ
// -----

// The next four-line record, or false at a clean end of file.
inline bool bio_next_fastq(BioLineReader& in, BioFastqEntry& e) {
    std::string h, s, p, q;
    do {
        if (!in.next(h)) return false;
    } while (h.empty());
    if (h[0] != '@') throw in.error("expected a FASTQ header starting with '@'");
    if (!in.next(s)) throw in.error("truncated FASTQ record");
    if (!in.next(p) || p.empty() || p[0] != '+') throw in.error("expected a '+' line");
    if (!in.next(q)) throw in.error("truncated FASTQ record");
    if (q.size() != s.size()) throw in.error("quality length differs from sequence length");
    e.header = h.substr(1);
    e.bases.bases.assign(s.begin(), s.end());
    e.bases.quals.assign(q.begin(), q.end());
    return true;
}

inline void bio_fastq_text(std::string& s, const BioFastqEntry& e) {
    if (e.bases.bases.size() != e.bases.quals.size()) {
        throw std::runtime_error("FASTQ entry '" + e.header + "': quality length differs from sequence length");
    }
    s += '@';
    s += e.header;
    s += '\n';
    s.append(reinterpret_cast<const char*>(e.bases.bases.data()), e.bases.bases.size());
    s += "\n+\n";
    s.append(reinterpret_cast<const char*>(e.bases.quals.data()), e.bases.quals.size());
    s += '\n';
}

inline std::string bio_format_fastq(const std::vector<BioFastqEntry>& batch) {
    std::string s;
    for (const auto& e : batch) bio_fastq_text(s, e);
    return s;
}

template <class Sink>
void bio_read_fastq(const std::string& path, Sink sink) {
    BioLineReader in(path);
    BioBatcher<BioFastqEntry, Sink> out(sink, BIO_BATCH_BYTES);
    BioFastqEntry e;
    while (bio_next_fastq(in, e)) {
        const size_t w = e.header.size() + 2 * e.bases.bases.size();
        out.add(std::move(e), w);
        e = BioFastqEntry{};
    }
    out.flush();
}

template <class Producer>
void bio_write_fastq(const std::string& path, Producer produce) {
    BioWriter out(path);
    bio_drain(produce, [&](const std::vector<BioFastqEntry>& batch) { out.write(bio_format_fastq(batch)); });
    out.close();
}

template <class Sink>
void bio_read_fastq_pairs(const std::string& path1, const std::string& path2, Sink sink) {
    using Pair = std::tuple<BioFastqEntry, BioFastqEntry>;
    BioLineReader in1(path1), in2(path2);
    BioBatcher<Pair, Sink> out(sink, BIO_BATCH_BYTES);
    while (true) {
        BioFastqEntry a, b;
        const bool more1 = bio_next_fastq(in1, a);
        const bool more2 = bio_next_fastq(in2, b);
        if (more1 != more2) {
            throw std::runtime_error("mate files " + path1 + " and " + path2 + " hold different numbers of reads");
        }
        if (!more1) break;
        const size_t w = a.header.size() + b.header.size() + 2 * (a.bases.bases.size() + b.bases.bases.size());
        out.add(Pair{std::move(a), std::move(b)}, w);
    }
    out.flush();
}

template <class Producer>
void bio_write_fastq_pairs(const std::string& path1, const std::string& path2, Producer produce) {
    BioWriter out1(path1), out2(path2);
    bio_drain(produce, [&](const std::vector<std::tuple<BioFastqEntry, BioFastqEntry>>& batch) {
        std::string s1, s2;
        for (const auto& [a, b] : batch) {
            bio_fastq_text(s1, a);
            bio_fastq_text(s2, b);
        }
        out1.write(s1);
        out2.write(s2);
    });
    out1.close();
    out2.close();
}


// GFF3 and GTF
// ------------

// One line of a GFF3 or GTF file.
struct BioGffLine {
    std::string seqid;
    std::string origin;
    std::string kind;
    BioSegment part;
    std::optional<double> score;
    BioPhase phase;
};

inline void bio_attr_add(BioAttrs& attrs, const std::string& key, std::vector<std::string> values) {
    for (auto& [k, vs] : attrs) {
        if (k == key) {
            vs.insert(vs.end(), values.begin(), values.end());
            return;
        }
    }
    attrs.emplace_back(key, std::move(values));
}

inline const std::vector<std::string>* bio_attr_get(const BioAttrs& attrs, const std::string& key) {
    for (const auto& [k, vs] : attrs) {
        if (k == key) return &vs;
    }
    return nullptr;
}

inline std::optional<std::string> bio_attr_one(const BioLineReader& in, const BioAttrs& attrs,
                                               const std::string& key) {
    const auto* vs = bio_attr_get(attrs, key);
    if (!vs || vs->empty()) return std::nullopt;
    if (vs->size() > 1) throw in.error("attribute " + key + " has several values");
    return (*vs)[0];
}

inline int bio_hex(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

// Decodes percent escapes; `fail` reports a malformed one.
template <class Fail>
std::string bio_percent_decode_with(const std::string& s, Fail fail) {
    std::string out;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] != '%') {
            out.push_back(s[i]);
            continue;
        }
        if (i + 2 >= s.size()) fail("truncated percent escape in '" + s + "'");
        const int hi = bio_hex(s[i + 1]);
        const int lo = bio_hex(s[i + 2]);
        if (hi < 0 || lo < 0) fail("bad percent escape in '" + s + "'");
        out.push_back(static_cast<char>(hi * 16 + lo));
        i += 2;
    }
    return out;
}

inline std::string bio_percent_decode(const BioLineReader& in, const std::string& s) {
    return bio_percent_decode_with(s, [&](const std::string& m) { throw in.error(m); });
}

// Escapes the characters GFF3 reserves in attribute keys and values.
inline std::string bio_percent_encode(const std::string& s) {
    static const char* hex = "0123456789ABCDEF";
    std::string out;
    for (unsigned char c : s) {
        if (c < 0x20 || c == 0x7f || c == '%' || c == ';' || c == '=' || c == '&' || c == ',') {
            out.push_back('%');
            out.push_back(hex[c >> 4]);
            out.push_back(hex[c & 15]);
        } else {
            out.push_back(static_cast<char>(c));
        }
    }
    return out;
}

// Parses a GFF3 attribute column (as `GffNode.attrs` holds it).
inline BioAttrs bio_attrs_parse(const std::string& col) {
    auto fail = [](const std::string& m) -> void { throw std::runtime_error("attributes: " + m); };
    BioAttrs attrs;
    if (col.empty() || col == ".") return attrs;
    for (const auto& field : bio_split(col, ';')) {
        const std::string f = bio_trim(field);
        if (f.empty()) continue;
        const size_t eq = f.find('=');
        if (eq == std::string::npos) fail("attribute without '=': '" + f + "'");
        std::vector<std::string> values;
        for (const auto& v : bio_split(f.substr(eq + 1), ',')) values.push_back(bio_percent_decode_with(v, fail));
        bio_attr_add(attrs, bio_percent_decode_with(f.substr(0, eq), fail), std::move(values));
    }
    return attrs;
}

// Attributes as a GFF3 attribute column.
inline std::string bio_attrs_text(const BioAttrs& attrs) {
    std::string s;
    for (const auto& [k, vs] : attrs) {
        if (!s.empty()) s += ";";
        s += bio_percent_encode(k) + "=";
        for (size_t i = 0; i < vs.size(); ++i) {
            if (i) s += ",";
            s += bio_percent_encode(vs[i]);
        }
    }
    return s;
}

inline std::vector<std::string> bio_attr_values(const std::string& key, const std::string& attrs) {
    const auto parsed = bio_attrs_parse(attrs);
    const auto* vs = bio_attr_get(parsed, key);
    return vs ? *vs : std::vector<std::string>{};
}

// Splits a GFF3 attribute column into its ID, its Parents and the rest of
// the column, which is kept as written. Escapes are checked everywhere.
inline std::string bio_split_gff3_attrs(const BioLineReader& in, const std::string& col, std::string& id,
                                        std::vector<std::string>& parents) {
    std::string rest;
    if (col == ".") return rest;
    bool haveId = false;
    for (const auto& field : bio_split(col, ';')) {
        const std::string f = bio_trim(field);
        if (f.empty()) continue;
        const size_t eq = f.find('=');
        if (eq == std::string::npos) throw in.error("attribute without '=': '" + f + "'");
        const std::string key = bio_percent_decode(in, f.substr(0, eq));
        std::vector<std::string> values;
        for (const auto& v : bio_split(f.substr(eq + 1), ',')) values.push_back(bio_percent_decode(in, v));
        if (key == "ID") {
            if (haveId || values.size() != 1) throw in.error("attribute ID has several values");
            id = values[0];
            haveId = true;
        } else if (key == "Parent") {
            parents.insert(parents.end(), values.begin(), values.end());
        } else {
            if (!rest.empty()) rest += ";";
            rest += f;
        }
    }
    return rest;
}

// GTF attributes: `key "value"; key value;`, a key repeated for several
// values.
inline BioAttrs bio_parse_gtf_attrs(const BioLineReader& in, const std::string& col) {
    BioAttrs attrs;
    size_t i = 0;
    const size_t n = col.size();
    auto skip_space = [&]() {
        while (i < n && (col[i] == ' ' || col[i] == '\t')) ++i;
    };
    while (true) {
        skip_space();
        while (i < n && col[i] == ';') {
            ++i;
            skip_space();
        }
        if (i >= n) break;
        const size_t k0 = i;
        while (i < n && col[i] != ' ' && col[i] != '\t' && col[i] != ';') ++i;
        const std::string key = col.substr(k0, i - k0);
        skip_space();
        std::string value;
        if (i < n && col[i] == '"') {
            const size_t close = col.find('"', i + 1);
            if (close == std::string::npos) throw in.error("unterminated quote in attributes");
            value = col.substr(i + 1, close - i - 1);
            i = close + 1;
        } else {
            const size_t v0 = i;
            while (i < n && col[i] != ';') ++i;
            value = bio_trim(col.substr(v0, i - v0));
        }
        if (key.empty()) throw in.error("attribute without a key");
        bio_attr_add(attrs, key, {value});
        skip_space();
        if (i < n && col[i] != ';') throw in.error("expected ';' after attribute " + key);
    }
    return attrs;
}

inline char bio_phase_char(BioPhase p) {
    switch (p) {
        case BioPhase::Phase0: return '0';
        case BioPhase::Phase1: return '1';
        case BioPhase::Phase2: return '2';
        default: return '.';
    }
}

// The nine columns of a feature line, attributes left unparsed.
inline BioGffLine bio_parse_gff_columns(const BioLineReader& in, const std::string& line, std::string& attrCol) {
    const auto cols = bio_split(line, '\t');
    if (cols.size() != 9) throw in.error("expected 9 tab-separated columns, found " + std::to_string(cols.size()));
    BioGffLine g;
    g.seqid = cols[0];
    g.origin = cols[1] == "." ? "" : cols[1];
    g.kind = cols[2];
    const uint64_t start = bio_parse_u64(in, cols[3], "start");
    const uint64_t end = bio_parse_u64(in, cols[4], "end");
    if (start == 0) throw in.error("start is 0; GFF coordinates are 1-based");
    if (start > end + 1) throw in.error("start is past end");
    g.part = BioSegment{start - 1, end, bio_parse_strand(in, cols[6])};
    g.score = bio_parse_score(in, cols[5]);
    if (cols[7] == ".") g.phase = BioPhase::PhaseNone;
    else if (cols[7] == "0") g.phase = BioPhase::Phase0;
    else if (cols[7] == "1") g.phase = BioPhase::Phase1;
    else if (cols[7] == "2") g.phase = BioPhase::Phase2;
    else throw in.error("bad phase '" + cols[7] + "'");
    attrCol = cols[8];
    return g;
}

// A node under construction: the node, one segment, score and phase per
// line, and where it sits.
struct BioGffBuild {
    std::string seqid;
    BioGffNode node;
    std::vector<BioSegment> parts;
    std::vector<std::string> parents;
    std::vector<size_t> children;
};

// The feature of a built node: parts in transcript order, with scores and
// phases in the same order, each dropped when every line's is '.'.
inline BioFeature<BioGffNode> bio_built_feature(const BioGffBuild& b) {
    const auto perm = bio_transcript_order(b.parts);
    BioFeature<BioGffNode> f{b.node, bio_permute(b.parts, perm), {}};
    f.node.scores = bio_permute(b.node.scores, perm);
    f.node.phases = bio_permute(b.node.phases, perm);
    if (std::all_of(f.node.scores.begin(), f.node.scores.end(), [](const auto& x) { return !x; })) f.node.scores.clear();
    if (std::all_of(f.node.phases.begin(), f.node.phases.end(), [](BioPhase x) { return x == BioPhase::PhaseNone; }))
        f.node.phases.clear();
    return f;
}

// Builds the trees of a set of nodes whose Parent references resolve among
// themselves, in order of first appearance, copying a node under each of its
// parents.
inline void bio_gff_trees(const std::string& path, std::vector<BioGffBuild>& nodes,
                          const std::unordered_map<std::string, size_t>& byId,
                          std::vector<BioGffRoot>& roots) {
    std::vector<size_t> tops;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].parents.empty()) {
            tops.push_back(i);
            continue;
        }
        for (const auto& p : nodes[i].parents) {
            const auto it = byId.find(p);
            if (it == byId.end()) {
                throw std::runtime_error(path + ": Parent '" + p + "' of a " + nodes[i].node.kind +
                                         " feature is not defined");
            }
            if (nodes[it->second].seqid != nodes[i].seqid) {
                throw std::runtime_error(path + ": feature tree '" + p + "' spans several sequences");
            }
            nodes[it->second].children.push_back(i);
        }
    }
    std::vector<int> state(nodes.size(), 0);  // 0 new, 1 on path, 2 done
    std::function<BioFeature<BioGffNode>(size_t)> build = [&](size_t i) {
        if (state[i] == 1) throw std::runtime_error(path + ": Parent references form a cycle");
        state[i] = 1;
        BioFeature<BioGffNode> f = bio_built_feature(nodes[i]);
        for (size_t c : nodes[i].children) f.children.push_back(build(c));
        state[i] = 2;
        return f;
    };
    for (size_t t : tops) roots.push_back(BioGffRoot{nodes[t].seqid, build(t)});
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (state[i] == 0) throw std::runtime_error(path + ": Parent references form a cycle");
    }
}

// Reads GFF3. Trees are complete at a "###" line or at the end of the file;
// a file without "###" lines is held in memory until its end.
template <class Sink>
void bio_read_gff3(const std::string& path, Sink sink) {
    BioLineReader in(path);
    BioBatcher<BioGffRoot, Sink> out(sink, BIO_BATCH_TREES);
    std::vector<BioGffBuild> nodes;
    std::unordered_map<std::string, size_t> byId;
    auto flush_trees = [&]() {
        std::vector<BioGffRoot> roots;
        bio_gff_trees(path, nodes, byId, roots);
        for (auto& r : roots) out.add(std::move(r), 1);
        nodes.clear();
        byId.clear();
    };
    std::string line, attrCol;
    while (in.next(line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (line == "###") flush_trees();
            if (line.rfind("##FASTA", 0) == 0) break;
            continue;
        }
        BioGffLine g = bio_parse_gff_columns(in, line, attrCol);
        std::string id;
        std::vector<std::string> parents;
        std::string attrs = bio_split_gff3_attrs(in, attrCol, id, parents);
        if (!id.empty()) {
            const auto it = byId.find(id);
            if (it != byId.end()) {
                // Another line of a discontinuous feature.
                BioGffBuild& b = nodes[it->second];
                if (b.seqid != g.seqid || b.node.kind != g.kind) {
                    throw in.error("ID '" + id + "' is reused by a different feature");
                }
                b.parts.push_back(g.part);
                b.node.scores.push_back(g.score);
                b.node.phases.push_back(g.phase);
                continue;
            }
            byId[id] = nodes.size();
        }
        BioGffBuild b;
        b.seqid = g.seqid;
        b.node = BioGffNode{g.kind, g.origin, id, std::move(attrs), {g.score}, {g.phase}, false};
        b.parts = {g.part};
        b.parents = std::move(parents);
        nodes.push_back(std::move(b));
    }
    flush_trees();
    out.flush();
}

// Checks that a node's scores and phases are empty or one per part.
inline void bio_check_parallel(const BioGffNode& n, size_t nparts) {
    if ((!n.scores.empty() && n.scores.size() != nparts) || (!n.phases.empty() && n.phases.size() != nparts)) {
        throw std::runtime_error("a " + n.kind + " node's scores and phases do not match its parts");
    }
}

inline std::string bio_node_score(const BioGffNode& n, size_t i) {
    return n.scores.empty() ? std::string(".") : bio_format_score(n.scores[i]);
}

inline char bio_node_phase(const BioGffNode& n, size_t i) {
    return n.phases.empty() ? '.' : bio_phase_char(n.phases[i]);
}

// One tree as GFF3 lines, followed by "###". A node copied under several
// parents is written once, listing all of them. A node without an ID that
// has children gets `prefix` followed by a count kept in `fresh`.
inline void bio_gff3_text(std::string& s, const BioGffRoot& root, const std::string& prefix, size_t& fresh) {
    std::map<std::string, std::vector<std::string>> parentsOf;
    std::vector<std::pair<const BioFeature<BioGffNode>*, std::string>> order;
    std::function<void(const BioFeature<BioGffNode>&, const std::optional<std::string>&)> visit =
        [&](const BioFeature<BioGffNode>& f, const std::optional<std::string>& parent) {
            std::string id = f.node.id;
            if (id.empty() && !f.children.empty()) id = prefix + std::to_string(++fresh);
            const bool seen = !id.empty() && parentsOf.count(id);
            if (!id.empty()) {
                auto& ps = parentsOf[id];
                if (parent) ps.push_back(*parent);
            }
            if (!seen) order.emplace_back(&f, id);
            if (id.empty() && parent) {
                // No ID: remember the single parent in the order entry.
                order.back().second = "\x01" + *parent;
            }
            if (!seen) {
                for (const auto& c : f.children) visit(c, id.empty() ? std::nullopt : std::optional<std::string>(id));
            }
        };
    visit(root.feature, std::nullopt);
    for (const auto& [f, key] : order) {
        const BioGffNode& n = f->node;
        bio_check_parallel(n, f->parts.size());
        std::string attrs;
        auto append = [&](const std::string& a) {
            if (a.empty()) return;
            if (!attrs.empty()) attrs += ";";
            attrs += a;
        };
        if (!key.empty() && key[0] == '\x01') {
            append("Parent=" + bio_percent_encode(key.substr(1)));
        } else if (!key.empty()) {
            append("ID=" + bio_percent_encode(key));
            const auto& ps = parentsOf[key];
            if (!ps.empty()) {
                std::string pa = "Parent=";
                for (size_t i = 0; i < ps.size(); ++i) {
                    if (i) pa += ",";
                    pa += bio_percent_encode(ps[i]);
                }
                append(pa);
            }
        }
        append(n.attrs);
        if (attrs.empty()) attrs = ".";
        const std::string origin = n.origin.empty() ? "." : n.origin;
        for (size_t i = 0; i < f->parts.size(); ++i) {
            const BioSegment& p = f->parts[i];
            s += root.seqid + "\t" + origin + "\t" + n.kind + "\t" + std::to_string(p.start + 1) + "\t" +
                 std::to_string(p.end) + "\t" + bio_node_score(n, i) + "\t" + bio_strand_char(p.strand) + "\t" +
                 bio_node_phase(n, i) + "\t" + attrs + "\n";
        }
    }
    s += "###\n";
}

// A batch of trees as GFF3 text. The header is written for the batch at
// offset 0; generated IDs carry the offset, so batches never collide.
inline std::string bio_format_gff3(int offset, const std::vector<BioGffRoot>& batch) {
    std::string s = offset == 0 ? "##gff-version 3\n" : "";
    const std::string prefix = "bio_node_" + std::to_string(offset) + "_";
    size_t fresh = 0;
    for (const auto& root : batch) bio_gff3_text(s, root, prefix, fresh);
    return s;
}

template <class Producer>
void bio_write_gff3(const std::string& path, Producer produce) {
    BioWriter out(path);
    out.write(std::string("##gff-version 3\n"));
    size_t fresh = 0;
    bio_drain(produce, [&](const std::vector<BioGffRoot>& batch) {
        std::string s;
        for (const auto& root : batch) bio_gff3_text(s, root, "bio_node_", fresh);
        out.write(s);
    });
    out.close();
}

// A gene assembled from GTF lines.
struct BioGtfGene {
    std::string seqid;
    std::optional<BioGffBuild> gene;
    std::vector<std::string> txOrder;
    std::unordered_map<std::string, BioGffBuild> tx;
    std::unordered_map<std::string, std::vector<BioGffBuild>> txKids;
    std::vector<BioGffBuild> geneKids;
};

// The span covering a set of parts: strand shared by all, else Unknown.
inline BioSegment bio_envelope(const std::vector<BioSegment>& parts) {
    uint64_t s = UINT64_MAX, e = 0;
    for (const auto& p : parts) {
        s = std::min(s, p.start);
        e = std::max(e, p.end);
    }
    BioStrand st = parts.empty() ? BioStrand::Unknown : parts[0].strand;
    for (const auto& p : parts) {
        if (p.strand != st) st = BioStrand::Unknown;
    }
    return BioSegment{parts.empty() ? 0 : s, e, st};
}

// A gene or transcript node a GTF implies but does not write, spanning
// `parts`.
inline BioFeature<BioGffNode> bio_created(const std::string& kind, const std::string& id, const BioAttrs& attrs,
                                          const std::string& origin, const std::vector<BioSegment>& parts) {
    BioGffNode n{kind, origin, id, bio_attrs_text(attrs), {}, {}, true};
    return bio_feature(std::move(n), {bio_envelope(parts)});
}

inline BioGffRoot bio_gtf_tree(const std::string& geneId, BioGtfGene& g) {
    std::vector<BioSegment> geneParts;
    std::string origin;
    std::vector<BioFeature<BioGffNode>> txs;
    for (const auto& t : g.txOrder) {
        const auto& kids = g.txKids[t];
        std::vector<BioSegment> txParts;
        for (const auto& k : kids) {
            txParts.insert(txParts.end(), k.parts.begin(), k.parts.end());
            if (origin.empty()) origin = k.node.origin;
        }
        BioFeature<BioGffNode> tf;
        const auto it = g.tx.find(t);
        if (it != g.tx.end()) {
            tf = bio_built_feature(it->second);
            if (origin.empty()) origin = it->second.node.origin;
        } else {
            BioAttrs a;
            a.emplace_back("gene_id", std::vector<std::string>{geneId});
            a.emplace_back("transcript_id", std::vector<std::string>{t});
            tf = bio_created("transcript", t, a, kids.empty() ? "" : kids[0].node.origin, txParts);
        }
        for (const auto& k : kids) tf.children.push_back(bio_built_feature(k));
        geneParts.insert(geneParts.end(), tf.parts.begin(), tf.parts.end());
        txs.push_back(std::move(tf));
    }
    for (const auto& k : g.geneKids) {
        geneParts.insert(geneParts.end(), k.parts.begin(), k.parts.end());
        if (origin.empty()) origin = k.node.origin;
    }
    BioFeature<BioGffNode> gf;
    if (g.gene) {
        gf = bio_built_feature(*g.gene);
    } else {
        BioAttrs a;
        a.emplace_back("gene_id", std::vector<std::string>{geneId});
        gf = bio_created("gene", geneId, a, origin, geneParts);
    }
    for (auto& t : txs) gf.children.push_back(std::move(t));
    for (const auto& k : g.geneKids) gf.children.push_back(bio_built_feature(k));
    return BioGffRoot{g.seqid, std::move(gf)};
}

// Reads GTF, keeping the lines `keep(seqid, attrs)` accepts. A sequence's
// genes are held until the sequence ends, so a GTF sorted by position, whose
// genes interleave, reads like one grouped by gene; a gene_id that reappears
// on a later sequence is an error.
// Builds a GTF's trees from its lines, which `next(line, attrs)` yields
// until it returns false. A sequence's genes are held until the sequence
// ends, so a GTF sorted by position, whose genes interleave, reads like one
// grouped by gene; a gene_id that reappears on a later sequence is an error.
// `err` places an error message.
template <class Next, class Err, class Sink>
void bio_gtf_build(Next next, Err err, Sink sink) {
    BioBatcher<BioGffRoot, Sink> out(sink, BIO_BATCH_TREES);
    std::unordered_set<std::string> done;
    std::vector<std::pair<std::string, BioGtfGene>> pending;
    std::unordered_map<std::string, size_t> slot;
    std::optional<std::string> curSeq;
    auto flush = [&]() {
        for (auto& [id, gene] : pending) {
            out.add(bio_gtf_tree(id, gene), 1);
            done.insert(id);
        }
        pending.clear();
        slot.clear();
    };
    BioGffLine g;
    BioAttrs attrs;
    auto one = [&](const std::string& key) -> std::optional<std::string> {
        const auto* vs = bio_attr_get(attrs, key);
        if (!vs || vs->empty()) return std::nullopt;
        if (vs->size() > 1) throw err("attribute " + key + " has several values");
        return (*vs)[0];
    };
    while (next(g, attrs)) {
        const std::optional<std::string> geneId = one("gene_id");
        const std::optional<std::string> txId = one("transcript_id");
        if (!geneId || geneId->empty()) throw err("GTF line without a gene_id");
        if (!curSeq || *curSeq != g.seqid) {
            if (slot.count(*geneId)) throw err("gene " + *geneId + " spans several sequences");
            flush();
            curSeq = g.seqid;
        }
        size_t k;
        auto it = slot.find(*geneId);
        if (it != slot.end()) {
            k = it->second;
        } else {
            if (done.count(*geneId)) throw err("lines of gene " + *geneId + " are not contiguous");
            BioGtfGene fresh;
            fresh.seqid = g.seqid;
            pending.emplace_back(*geneId, std::move(fresh));
            k = pending.size() - 1;
            slot[*geneId] = k;
        }
        BioGtfGene& cur = pending[k].second;
        BioGffBuild b;
        b.seqid = g.seqid;
        b.parts = {g.part};
        const bool isGene = g.kind == "gene";
        const bool isTx = g.kind == "transcript";
        std::string id;
        if (isGene) id = *geneId;
        if (isTx && txId) id = *txId;
        b.node = BioGffNode{g.kind, g.origin, id, bio_attrs_text(attrs), {g.score}, {g.phase}, false};
        if (isGene) {
            if (cur.gene) throw err("gene " + *geneId + " has two gene lines");
            cur.gene = std::move(b);
        } else if (isTx && txId && !txId->empty()) {
            if (cur.tx.count(*txId)) throw err("transcript " + *txId + " has two transcript lines");
            if (!cur.txKids.count(*txId)) cur.txOrder.push_back(*txId);
            cur.txKids[*txId];
            cur.tx[*txId] = std::move(b);
        } else if (txId && !txId->empty()) {
            if (!cur.txKids.count(*txId)) cur.txOrder.push_back(*txId);
            cur.txKids[*txId].push_back(std::move(b));
        } else {
            cur.geneKids.push_back(std::move(b));
        }
    }
    flush();
    out.flush();
}


// Reads GTF, keeping the lines `keep(seqid, attrs)` accepts.
template <class Keep, class Sink>
void bio_read_gtf_with(const std::string& path, Keep keep, Sink sink) {
    BioLineReader in(path);
    std::string line, attrCol;
    bio_gtf_build(
        [&](BioGffLine& g, BioAttrs& attrs) {
            while (in.next(line)) {
                if (line.empty() || line[0] == '#') continue;
                g = bio_parse_gff_columns(in, line, attrCol);
                attrs = bio_parse_gtf_attrs(in, attrCol);
                if (keep(g.seqid, attrs)) return true;
            }
            return false;
        },
        [&](const std::string& m) { return in.error(m); },
        sink);
}

template <class Sink>
void bio_read_gtf(const std::string& path, Sink sink) {
    bio_read_gtf_with(path, [](const std::string&, const BioAttrs&) { return true; }, sink);
}

// A GTF stream keeping only lines on the named sequences and, when `checkTx`
// holds, with a non-empty transcript_id (nf-core's gtffilter).
template <class Sink>
void bio_read_gtf_filtered(const std::vector<std::string>& names, bool checkTx, const std::string& path, Sink sink) {
    const std::unordered_set<std::string> keep(names.begin(), names.end());
    bio_read_gtf_with(
        path,
        [&](const std::string& seqid, const BioAttrs& attrs) {
            if (!keep.count(seqid)) return false;
            if (!checkTx) return true;
            const std::vector<std::string>* v = bio_attr_get(attrs, "transcript_id");
            return v && std::any_of(v->begin(), v->end(), [](const std::string& t) { return !t.empty(); });
        },
        sink);
}

struct BioGtfLine {
    std::string seqid;
    std::string origin;
    std::string kind;
    BioSegment part;
    std::optional<double> score;
    BioPhase phase;
    std::string attrs;
};

template <class Sink>
void bio_read_gtf_lines(const std::string& path, Sink sink) {
    BioLineReader in(path);
    BioBatcher<BioGtfLine, Sink> out(sink, BIO_BATCH_TREES);
    std::string line, attrCol;
    while (in.next(line)) {
        if (line.empty() || line[0] == '#') continue;
        BioGffLine g = bio_parse_gff_columns(in, line, attrCol);
        BioAttrs attrs = bio_parse_gtf_attrs(in, attrCol);
        out.add(BioGtfLine{g.seqid, g.origin, g.kind, g.part, g.score, g.phase, bio_attrs_text(attrs)}, 1);
    }
    out.flush();
}

inline std::vector<BioGffRoot> bio_gtf_trees(const std::vector<BioGtfLine>& lines) {
    std::vector<BioGffRoot> roots;
    size_t i = 0;
    bio_gtf_build(
        [&](BioGffLine& g, BioAttrs& attrs) {
            if (i >= lines.size()) return false;
            const BioGtfLine& l = lines[i++];
            g = BioGffLine{l.seqid, l.origin, l.kind, l.part, l.score, l.phase};
            attrs = bio_attrs_parse(l.attrs);
            return true;
        },
        [&](const std::string& m) { return std::runtime_error("GTF line " + std::to_string(i) + ": " + m); },
        [&](const std::vector<BioGffRoot>& xs) { roots.insert(roots.end(), xs.begin(), xs.end()); });
    return roots;
}

inline std::vector<BioGtfLine> bio_filter_gtf_lines(const std::vector<std::string>& names, bool checkTx,
                                                     const std::vector<BioGtfLine>& lines) {
    const std::unordered_set<std::string> keep(names.begin(), names.end());
    std::vector<BioGtfLine> out;
    for (const auto& l : lines) {
        if (!keep.count(l.seqid)) continue;
        if (checkTx) {
            const BioAttrs attrs = bio_attrs_parse(l.attrs);
            const std::vector<std::string>* v = bio_attr_get(attrs, "transcript_id");
            if (!v || !std::any_of(v->begin(), v->end(), [](const std::string& t) { return !t.empty(); })) continue;
        }
        out.push_back(l);
    }
    return out;
}

inline std::string bio_gtf_attrs(const BioAttrs& attrs) {
    std::string s;
    for (const auto& [k, vs] : attrs) {
        for (const auto& v : vs) {
            if (!s.empty()) s += " ";
            s += k + " \"" + v + "\";";
        }
    }
    return s;
}

// Writes GTF. Nodes the reader created are not written. GTF's hierarchy is
// gene, transcript, feature, so ids follow tree depth: a line's gene_id is
// its own gene_id attribute or its root's, and its transcript_id is its own
// transcript_id attribute or that of the root's child it lies under (for
// either, the node's id stands in for a missing attribute).
inline void bio_gtf_text(std::string& s, const BioGffRoot& root) {
    std::function<void(const BioFeature<BioGffNode>&, size_t, const std::optional<std::string>&,
                       const std::optional<std::string>&)>
        visit = [&](const BioFeature<BioGffNode>& f, size_t depth, const std::optional<std::string>& geneId,
                    const std::optional<std::string>& txId) {
            const BioGffNode& n = f.node;
            const BioAttrs nattrs = bio_attrs_parse(n.attrs);
            const std::optional<std::string> nid = n.id.empty() ? std::nullopt : std::optional<std::string>(n.id);
            std::optional<std::string> g = geneId, t = txId;
            if (const auto* v = bio_attr_get(nattrs, "gene_id"); v && !v->empty()) g = (*v)[0];
            else if (depth == 0) g = nid;
            if (const auto* v = bio_attr_get(nattrs, "transcript_id"); v && !v->empty()) t = (*v)[0];
            else if (depth == 1) t = nid;
            if (!n.created) {
                if (!g) throw std::runtime_error("cannot write a " + n.kind + " feature to GTF without a gene_id");
                bio_check_parallel(n, f.parts.size());
                BioAttrs a;
                a.emplace_back("gene_id", std::vector<std::string>{*g});
                if (t && depth > 0) a.emplace_back("transcript_id", std::vector<std::string>{*t});
                for (const auto& kv : nattrs) {
                    const auto& k = std::get<0>(kv);
                    if (k != "gene_id" && k != "transcript_id") a.push_back(kv);
                }
                const std::string attrs = bio_gtf_attrs(a);
                const std::string origin = n.origin.empty() ? "." : n.origin;
                for (size_t i = 0; i < f.parts.size(); ++i) {
                    const BioSegment& p = f.parts[i];
                    s += root.seqid + "\t" + origin + "\t" + n.kind + "\t" + std::to_string(p.start + 1) + "\t" +
                         std::to_string(p.end) + "\t" + bio_node_score(n, i) + "\t" + bio_strand_char(p.strand) +
                         "\t" + bio_node_phase(n, i) + "\t" + attrs + "\n";
                }
            }
            for (const auto& c : f.children) visit(c, depth + 1, g, t);
        };
    visit(root.feature, 0, std::nullopt, std::nullopt);
}

inline std::string bio_format_gtf(const std::vector<BioGffRoot>& batch) {
    std::string s;
    for (const auto& root : batch) bio_gtf_text(s, root);
    return s;
}

template <class Producer>
void bio_write_gtf(const std::string& path, Producer produce) {
    BioWriter out(path);
    bio_drain(produce, [&](const std::vector<BioGffRoot>& batch) { out.write(bio_format_gtf(batch)); });
    out.close();
}


// Gene models
// -----------

inline bool bio_is_utr(const std::string& k) { return k == "five_prime_UTR" || k == "three_prime_UTR" || k == "UTR"; }
inline bool bio_is_codon(const std::string& k) { return k == "start_codon" || k == "stop_codon"; }

inline std::vector<BioSpan> bio_sorted_spans(std::vector<BioSpan> xs) {
    std::sort(xs.begin(), xs.end());
    return xs;
}

// The spans of `n` bases of `spans` (given in transcript order) starting
// `skip` bases in, in transcript order; fewer when they run out.
inline std::vector<BioSpan> bio_take_bases(const std::vector<BioSpan>& spans, bool minus, uint64_t skip, uint64_t n) {
    std::vector<BioSpan> out;
    for (const auto& [s, e] : spans) {
        if (n == 0) break;
        const uint64_t len = e - s;
        if (skip >= len) {
            skip -= len;
            continue;
        }
        const uint64_t take = std::min(n, len - skip);
        if (minus) out.emplace_back(e - skip - take, e - skip);
        else out.emplace_back(s + skip, s + skip + take);
        n -= take;
        skip = 0;
    }
    return out;
}

// The UTR spans of a coding transcript on either side of its CDS, each in
// transcript order.
inline void bio_utr_spans(const std::vector<BioSpan>& exonsAsc, BioSpan cds, bool minus, std::vector<BioSpan>& five,
                          std::vector<BioSpan>& three) {
    const auto [cs, ce] = cds;
    std::vector<BioSpan> lo, hi;
    for (const auto& [s, e] : exonsAsc) {
        if (s < cs) lo.emplace_back(s, std::min(e, cs));
        if (e > ce) hi.emplace_back(std::max(s, ce), e);
    }
    if (minus) {
        std::reverse(lo.begin(), lo.end());
        std::reverse(hi.begin(), hi.end());
        five = hi;
        three = lo;
    } else {
        five = lo;
        three = hi;
    }
}

// Phases of CDS parts in transcript order, derived from the first.
inline std::vector<BioPhase> bio_derived_phases(const std::vector<BioSpan>& partsTx, BioPhase first) {
    std::vector<BioPhase> out;
    if (first == BioPhase::PhaseNone) return out;
    uint64_t done = 0;
    const uint64_t p0 = static_cast<uint64_t>(first);
    for (size_t i = 0; i < partsTx.size(); ++i) {
        if (i == 0) {
            out.push_back(first);
        } else {
            const uint64_t left = (done + 3 - p0 % 3) % 3;  // (done - p0) mod 3
            out.push_back(static_cast<BioPhase>((3 - left) % 3));
        }
        done += std::get<1>(partsTx[i]) - std::get<0>(partsTx[i]);
    }
    return out;
}

inline bool bio_any_flag(const BioComplexity& c) {
    return c.mixedStrand || c.splitFeature || c.missingExons || c.overlappingExons || c.spanMismatch ||
           c.multipleCds || c.cdsNotInExons || c.phaseShift || c.underivedLines || c.extraChildren ||
           c.extraColumns || c.anonymousTranscript;
}

template <class N>
bool bio_tree_has(const BioFeature<N>& f, const std::string& kind) {
    if (f.node.kind == kind) return true;
    for (const auto& c : f.children)
        if (bio_tree_has(c, kind)) return true;
    return false;
}

// Collects every segment of a tree.
inline void bio_tree_segments(const BioFeature<BioGffNode>& f, std::vector<BioSegment>& out) {
    out.insert(out.end(), f.parts.begin(), f.parts.end());
    for (const auto& c : f.children) bio_tree_segments(c, out);
}

// Checks one transcript -- `t`'s children, with `t` itself absent for a gene
// whose CDS lines sit directly under it -- and builds its model.
inline BioGeneTx bio_check_tx(const BioFeature<BioGffNode>* t, const std::vector<const BioFeature<BioGffNode>*>& kids,
                              const BioGffNode& gene, bool minus, BioComplexity& c) {
    BioGeneTx tx;
    if (t) {
        tx.id = t->node.id;
        tx.kind = t->node.kind;
        tx.origin = t->node.origin;
        tx.attrs = t->node.attrs;
        if (tx.id.empty()) c.anonymousTranscript = true;
        if (t->parts.size() != 1) c.splitFeature = true;
    } else {
        tx.origin = kids.empty() ? gene.origin : kids[0]->node.origin;
    }
    std::vector<const BioFeature<BioGffNode>*> cdsNodes;
    std::vector<BioSpan> five, three, generic, start, stop;
    for (const auto* k : kids) {
        const BioGffNode& n = k->node;
        if (!k->children.empty()) c.extraChildren = true;
        if (!n.scores.empty()) c.extraColumns = true;
        if (n.origin != tx.origin) c.extraColumns = true;
        if (!n.phases.empty() && n.kind != "CDS") {
            // A codon split across an intron has phases like a CDS's, from 0.
            if (!bio_is_codon(n.kind) || bio_derived_phases(bio_spans(*k), BioPhase::Phase0) != n.phases)
                c.extraColumns = true;
        }
        if (n.kind == "exon") {
            if (k->parts.size() != 1) c.splitFeature = true;
            const auto sp = bio_spans(*k);
            tx.exons.insert(tx.exons.end(), sp.begin(), sp.end());
        } else if (n.kind == "CDS") {
            cdsNodes.push_back(k);
        } else if (n.kind == "five_prime_UTR") {
            const auto sp = bio_spans(*k);
            five.insert(five.end(), sp.begin(), sp.end());
        } else if (n.kind == "three_prime_UTR") {
            const auto sp = bio_spans(*k);
            three.insert(three.end(), sp.begin(), sp.end());
        } else if (n.kind == "UTR") {
            const auto sp = bio_spans(*k);
            generic.insert(generic.end(), sp.begin(), sp.end());
        } else if (n.kind == "start_codon") {
            const auto sp = bio_spans(*k);
            start.insert(start.end(), sp.begin(), sp.end());
        } else if (n.kind == "stop_codon") {
            const auto sp = bio_spans(*k);
            stop.insert(stop.end(), sp.begin(), sp.end());
        } else {
            c.extraChildren = true;
        }
    }
    if (!t) {
        // No transcript line: the exons are the CDS lines.
        for (const auto* k : cdsNodes) {
            const auto sp = bio_spans(*k);
            tx.exons.insert(tx.exons.end(), sp.begin(), sp.end());
        }
    }
    tx.exons = bio_sorted_spans(tx.exons);
    if (tx.exons.empty()) c.missingExons = true;
    for (size_t i = 1; i < tx.exons.size(); ++i)
        if (std::get<0>(tx.exons[i]) < std::get<1>(tx.exons[i - 1])) c.overlappingExons = true;
    if (t && t->parts.size() == 1 && !tx.exons.empty() &&
        BioSpan{t->parts[0].start, t->parts[0].end} != BioSpan{std::get<0>(tx.exons.front()), std::get<1>(tx.exons.back())})
        c.spanMismatch = true;
    // CDS lines without an ID are one CDS written line by line.
    BioFeature<BioGffNode> mergedCds;
    if (cdsNodes.size() > 1 &&
        std::all_of(cdsNodes.begin(), cdsNodes.end(), [](const auto* k) { return k->node.id.empty(); })) {
        BioGffBuild b;
        b.node = cdsNodes[0]->node;
        b.node.scores.clear();
        b.node.phases.clear();
        for (const auto* k : cdsNodes) {
            for (size_t i = 0; i < k->parts.size(); ++i) {
                b.parts.push_back(k->parts[i]);
                b.node.scores.push_back(k->node.scores.empty() ? std::nullopt : k->node.scores[i]);
                b.node.phases.push_back(k->node.phases.empty() ? BioPhase::PhaseNone : k->node.phases[i]);
            }
        }
        mergedCds = bio_built_feature(b);
        cdsNodes = {&mergedCds};
    }
    if (cdsNodes.size() > 1) c.multipleCds = true;
    std::vector<BioSpan> exonsTx = tx.exons;
    if (minus) std::reverse(exonsTx.begin(), exonsTx.end());
    if (cdsNodes.size() == 1 && !tx.exons.empty()) {
        const auto& cn = *cdsNodes[0];
        const auto cnSpans = bio_spans(cn);
        const auto cdsAsc = bio_sorted_spans(cnSpans);
        const BioSpan range{std::get<0>(cdsAsc.front()), std::get<1>(cdsAsc.back())};
        BioGeneTx probe;
        probe.exons = tx.exons;
        probe.cds = BioCds{cn.node.id, range, BioPhase::PhaseNone};
        if (bio_cds_spans(probe) != cdsAsc) c.cdsNotInExons = true;
        const BioPhase first = cn.node.phases.empty() ? BioPhase::PhaseNone : cn.node.phases[0];
        if (!cn.node.phases.empty() && bio_derived_phases(cnSpans, first) != cn.node.phases) c.phaseShift = true;
        tx.cds = BioCds{cn.node.id, range, first};
        // UTR and codon lines are kept only as what the model implies.
        std::vector<BioSpan> f5, f3;
        bio_utr_spans(tx.exons, range, minus, f5, f3);
        std::vector<BioSpan> both = f5;
        both.insert(both.end(), f3.begin(), f3.end());
        if ((!five.empty() || !three.empty()) &&
            (bio_sorted_spans(five) != bio_sorted_spans(f5) || bio_sorted_spans(three) != bio_sorted_spans(f3)))
            c.underivedLines = true;
        if (!generic.empty() && bio_sorted_spans(generic) != bio_sorted_spans(both)) c.underivedLines = true;
        const std::vector<BioSpan>& cdsTx = cnSpans;  // transcript order
        if (!start.empty() && bio_sorted_spans(start) != bio_sorted_spans(bio_take_bases(cdsTx, minus, 0, 3)))
            c.underivedLines = true;
        if (!stop.empty()) {
            uint64_t cdsLen = 0;
            for (const auto& [s, e] : cdsTx) cdsLen += e - s;
            const auto st = bio_sorted_spans(stop);
            const bool last3 = cdsLen >= 3 && st == bio_sorted_spans(bio_take_bases(cdsTx, minus, cdsLen - 3, 3));
            uint64_t before = 0;  // exon bases 5' of the CDS end
            for (const auto& [s, e] : exonsTx) {
                if (minus) before += e <= std::get<0>(range) ? 0 : e - std::max(s, std::get<0>(range));
                else before += s >= std::get<1>(range) ? 0 : std::min(e, std::get<1>(range)) - s;
            }
            const bool after3 = st == bio_sorted_spans(bio_take_bases(exonsTx, minus, before, 3));
            if (!last3 && !after3) c.underivedLines = true;
        }
    } else if (!five.empty() || !three.empty() || !generic.empty() || !start.empty() || !stop.empty()) {
        c.underivedLines = true;
    }
    return tx;
}

// Classifies one tree.
inline BioGffEntry bio_gff_entry(const BioGffRoot& r) {
    const auto& g = r.feature;
    auto other = [&]() { return BioGffEntry{std::make_shared<BioGffEntry_OtherFeature>(BioGffEntry_OtherFeature{r})}; };
    if (g.children.empty()) return other();
    std::vector<const BioFeature<BioGffNode>*> txs, leaves;
    for (const auto& ch : g.children) (ch.children.empty() ? leaves : txs).push_back(&ch);
    const bool implicit = txs.empty();
    if (implicit && !std::all_of(leaves.begin(), leaves.end(), [](const auto* k) { return k->node.kind == "CDS"; }))
        return other();
    BioComplexity c{};
    std::vector<BioSegment> all;
    bio_tree_segments(g, all);
    const BioStrand st = all[0].strand;
    for (const auto& p : all)
        if (p.strand != st) c.mixedStrand = true;
    const bool minus = st == BioStrand::Minus;
    if (g.parts.size() != 1) c.splitFeature = true;
    if (!g.node.scores.empty() || !g.node.phases.empty()) c.extraColumns = true;
    if (!implicit && !leaves.empty()) c.extraChildren = true;
    for (const auto* t : txs)
        if (!t->node.scores.empty() || !t->node.phases.empty()) c.extraColumns = true;
    BioGene gene;
    gene.seqid = r.seqid;
    gene.id = g.node.id;
    gene.kind = g.node.kind;
    gene.origin = g.node.origin;
    gene.attrs = g.node.attrs;
    gene.span = g.parts.empty() ? BioSpan{0, 0} : BioSpan{g.parts[0].start, g.parts[0].end};
    gene.strand = st;
    if (implicit) {
        gene.transcripts.push_back(bio_check_tx(nullptr, leaves, g.node, minus, c));
    } else {
        for (const auto* t : txs) {
            std::vector<const BioFeature<BioGffNode>*> kids;
            for (const auto& k : t->children) kids.push_back(&k);
            gene.transcripts.push_back(bio_check_tx(t, kids, g.node, minus, c));
        }
    }
    const bool coding = bio_tree_has(g, "CDS");
    if (bio_any_flag(c)) {
        if (coding) return BioGffEntry{std::make_shared<BioGffEntry_ComplexCodingGene>(BioGffEntry_ComplexCodingGene{c, r})};
        return BioGffEntry{std::make_shared<BioGffEntry_ComplexNonCodingGene>(BioGffEntry_ComplexNonCodingGene{c, r})};
    }
    if (coding) return BioGffEntry{std::make_shared<BioGffEntry_SimpleCodingGene>(BioGffEntry_SimpleCodingGene{std::move(gene)})};
    return BioGffEntry{std::make_shared<BioGffEntry_SimpleNonCodingGene>(BioGffEntry_SimpleNonCodingGene{std::move(gene)})};
}

inline std::vector<BioGffEntry> bio_gff_entries(const std::vector<BioGffRoot>& roots) {
    std::vector<BioGffEntry> out;
    out.reserve(roots.size());
    for (const auto& r : roots) out.push_back(bio_gff_entry(r));
    return out;
}

// A leaf node of a standard gene: no ID or attributes.
inline BioFeature<BioGffNode> bio_model_leaf(const std::string& kind, const std::string& id, const std::string& origin,
                                             BioStrand st, std::vector<BioSpan> partsTx,
                                             std::vector<BioPhase> phases = {}) {
    std::vector<BioSegment> segs;
    for (const auto& [a, b] : partsTx) segs.push_back(BioSegment{a, b, st});
    return BioFeature<BioGffNode>{BioGffNode{kind, origin, id, "", {}, std::move(phases), false}, std::move(segs), {}};
}

// The tree of a standard gene: the lines a GFF3 has for it.
inline BioGffRoot bio_gene_root(const BioGene& g) {
    const bool minus = g.strand == BioStrand::Minus;
    auto tx_order = [&](std::vector<BioSpan> xs) {
        if (minus) std::reverse(xs.begin(), xs.end());
        return xs;
    };
    BioFeature<BioGffNode> gf{BioGffNode{g.kind, g.origin, g.id, g.attrs, {}, {}, false},
                              {BioSegment{std::get<0>(g.span), std::get<1>(g.span), g.strand}}, {}};
    for (const auto& t : g.transcripts) {
        std::vector<BioFeature<BioGffNode>> kids;
        const auto exonsTx = tx_order(t.exons);
        if (!t.id.empty()) {
            for (const auto& e : exonsTx) kids.push_back(bio_model_leaf("exon", "", t.origin, g.strand, {e}));
        }
        if (t.cds) {
            const auto cdsTx = tx_order(bio_cds_spans(t));
            kids.push_back(bio_model_leaf("CDS", t.cds->id, t.origin, g.strand, cdsTx,
                                          bio_derived_phases(cdsTx, t.cds->phase)));
            if (!t.id.empty()) {
                std::vector<BioSpan> f5, f3;
                bio_utr_spans(t.exons, t.cds->range, minus, f5, f3);
                for (const auto& u : f5) kids.push_back(bio_model_leaf("five_prime_UTR", "", t.origin, g.strand, {u}));
                for (const auto& u : f3) kids.push_back(bio_model_leaf("three_prime_UTR", "", t.origin, g.strand, {u}));
            }
        }
        if (t.id.empty()) {
            for (auto& k : kids) gf.children.push_back(std::move(k));
            continue;
        }
        BioFeature<BioGffNode> tf{BioGffNode{t.kind, t.origin, t.id, t.attrs, {}, {}, false},
                                  {BioSegment{std::get<0>(t.exons.front()), std::get<1>(t.exons.back()), g.strand}},
                                  std::move(kids)};
        gf.children.push_back(std::move(tf));
    }
    return BioGffRoot{g.seqid, std::move(gf)};
}

inline std::vector<BioGffRoot> bio_entry_roots(const std::vector<BioGffEntry>& entries) {
    std::vector<BioGffRoot> out;
    out.reserve(entries.size());
    for (const auto& e : entries)
        bio_entry_visit(e, [&](const BioGene& g) { out.push_back(bio_gene_root(g)); },
                        [&](const BioGffRoot& r) { out.push_back(r); });
    return out;
}


// BED
// ---

// A score column: '.' is NaN.
inline double bio_parse_bed_score(const BioLineReader& in, const std::string& s) {
    if (s == ".") return std::numeric_limits<double>::quiet_NaN();
    return *bio_parse_score(in, s);
}

inline std::string bio_format_bed_score(double v) {
    return std::isnan(v) ? std::string(".") : bio_format_score(v);
}

// itemRgb: "0", or "r,g,b" with each in 0..255.
inline std::tuple<bool, uint8_t, uint8_t, uint8_t> bio_parse_rgb(const BioLineReader& in, const std::string& s) {
    if (s == "0") return {false, 0, 0, 0};
    const auto xs = bio_split(s, ',');
    if (xs.size() != 3) throw in.error("bad itemRgb '" + s + "'");
    uint8_t c[3];
    for (int i = 0; i < 3; ++i) {
        const uint64_t v = bio_parse_u64(in, xs[i], "itemRgb");
        if (v > 255) throw in.error("bad itemRgb '" + s + "'");
        c[i] = static_cast<uint8_t>(v);
    }
    return {true, c[0], c[1], c[2]};
}

inline std::string bio_format_rgb(const std::tuple<bool, uint8_t, uint8_t, uint8_t>& c) {
    if (!std::get<0>(c)) return "0";
    return std::to_string(std::get<1>(c)) + "," + std::to_string(std::get<2>(c)) + "," + std::to_string(std::get<3>(c));
}

// The first `n` tab-separated columns of `line`; fails on fewer.
inline void bio_bed_columns(const BioLineReader& in, const std::string& line, size_t n, std::vector<std::string>& cols) {
    cols.clear();
    size_t a = 0;
    while (cols.size() < n) {
        const size_t b = line.find('\t', a);
        cols.push_back(line.substr(a, b == std::string::npos ? std::string::npos : b - a));
        if (b == std::string::npos) break;
        a = b + 1;
    }
    if (cols.size() < n) throw in.error("expected at least " + std::to_string(n) + " tab-separated columns");
}

inline bool bio_bed_skip(const std::string& line) {
    return line.empty() || line[0] == '#' || line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0;
}

inline BioBed3 bio_bed3_of(const BioLineReader& in, const std::vector<std::string>& cols) {
    BioBed3 r{cols[0], bio_parse_u64(in, cols[1], "start"), bio_parse_u64(in, cols[2], "end")};
    if (r.start > r.end) throw in.error("start is past end");
    return r;
}

inline BioBed6 bio_bed6_of(const BioLineReader& in, const std::vector<std::string>& cols) {
    BioBed3 b = bio_bed3_of(in, cols);
    return BioBed6{std::move(b.chrom), b.start, b.end, cols[3], bio_parse_bed_score(in, cols[4]),
                   bio_parse_strand(in, cols[5])};
}

inline BioBed12 bio_bed12_of(const BioLineReader& in, const std::vector<std::string>& cols) {
    BioBed6 b = bio_bed6_of(in, cols);
    BioBed12 r{std::move(b.chrom), b.start, b.end, std::move(b.name), b.score, b.strand,
               bio_parse_u64(in, cols[6], "thickStart"), bio_parse_u64(in, cols[7], "thickEnd"),
               bio_parse_rgb(in, cols[8]), {}};
    if (r.thickStart > r.thickEnd || r.thickStart < r.start || r.thickEnd > r.end)
        throw in.error("thickStart and thickEnd must lie in order within the record");
    const uint64_t count = bio_parse_u64(in, cols[9], "blockCount");
    auto list = [&](const std::string& s, const char* what) {
        std::vector<uint64_t> v;
        for (const auto& x : bio_split(s, ',')) {
            if (!x.empty()) v.push_back(bio_parse_u64(in, x, what));
        }
        return v;
    };
    const auto sizes = list(cols[10], "blockSize");
    const auto starts = list(cols[11], "blockStart");
    if (count == 0 || sizes.size() != count || starts.size() != count)
        throw in.error("blockCount does not match blockSizes and blockStarts");
    uint64_t prev = r.start;
    for (uint64_t i = 0; i < count; ++i) {
        const uint64_t bs = r.start + starts[i], be = bs + sizes[i];
        if (i == 0 ? bs != r.start : bs < prev) throw in.error("blocks must start at chromStart, ascending and apart");
        r.blocks.emplace_back(bs, be);
        prev = be;
    }
    if (prev != r.end) throw in.error("the last block must end at chromEnd");
    return r;
}

template <class T, size_t N, class Of, class Sink>
void bio_read_bed_n(const std::string& path, Of of, Sink sink) {
    BioLineReader in(path);
    BioBatcher<T, Sink> out(sink, BIO_BATCH_TREES);
    std::string line;
    std::vector<std::string> cols;
    while (in.next(line)) {
        if (bio_bed_skip(line)) continue;
        bio_bed_columns(in, line, N, cols);
        out.add(of(in, cols), 1);
    }
    out.flush();
}

template <class Sink>
void bio_read_bed3(const std::string& path, Sink sink) {
    bio_read_bed_n<BioBed3, 3>(path, bio_bed3_of, sink);
}

template <class Sink>
void bio_read_bed6(const std::string& path, Sink sink) {
    bio_read_bed_n<BioBed6, 6>(path, bio_bed6_of, sink);
}

template <class Sink>
void bio_read_bed12(const std::string& path, Sink sink) {
    bio_read_bed_n<BioBed12, 12>(path, bio_bed12_of, sink);
}

inline void bio_bed3_cols(std::string& s, const std::string& chrom, uint64_t start, uint64_t end) {
    s += chrom;
    s += '\t';
    s += std::to_string(start);
    s += '\t';
    s += std::to_string(end);
}

inline void bio_bed6_cols(std::string& s, const std::string& chrom, uint64_t start, uint64_t end,
                          const std::string& name, double score, BioStrand strand) {
    bio_bed3_cols(s, chrom, start, end);
    s += '\t';
    s += name;
    s += '\t';
    s += bio_format_bed_score(score);
    s += '\t';
    s += bio_strand_char(strand);
}

inline std::string bio_format_bed3(const std::vector<BioBed3>& batch) {
    std::string s;
    for (const auto& r : batch) {
        bio_bed3_cols(s, r.chrom, r.start, r.end);
        s += '\n';
    }
    return s;
}

inline std::string bio_format_bed6(const std::vector<BioBed6>& batch) {
    std::string s;
    for (const auto& r : batch) {
        bio_bed6_cols(s, r.chrom, r.start, r.end, r.name, r.score, r.strand);
        s += '\n';
    }
    return s;
}

inline std::string bio_format_bed12(const std::vector<BioBed12>& batch) {
    std::string s;
    for (const auto& r : batch) {
        if (r.blocks.empty()) throw std::runtime_error("a BED12 record needs at least one block");
        bio_bed6_cols(s, r.chrom, r.start, r.end, r.name, r.score, r.strand);
        s += '\t' + std::to_string(r.thickStart) + '\t' + std::to_string(r.thickEnd) + '\t' + bio_format_rgb(r.rgb);
        std::string sizes, starts;
        for (const auto& [a, b] : r.blocks) {
            sizes += std::to_string(b - a) + ",";
            starts += std::to_string(a - r.start) + ",";
        }
        s += '\t' + std::to_string(r.blocks.size()) + '\t' + sizes + '\t' + starts + '\n';
    }
    return s;
}

template <class T, class Format, class Producer>
void bio_write_bed_n(const std::string& path, Format format, Producer produce) {
    BioWriter out(path);
    bio_drain(produce, [&](const std::vector<T>& batch) { out.write(format(batch)); });
    out.close();
}

template <class Producer>
void bio_write_bed3(const std::string& path, Producer produce) {
    bio_write_bed_n<BioBed3>(path, bio_format_bed3, produce);
}

template <class Producer>
void bio_write_bed6(const std::string& path, Producer produce) {
    bio_write_bed_n<BioBed6>(path, bio_format_bed6, produce);
}

template <class Producer>
void bio_write_bed12(const std::string& path, Producer produce) {
    bio_write_bed_n<BioBed12>(path, bio_format_bed12, produce);
}


// Reference products
// ------------------

// The first word of each header in a FASTA file, in order.
inline std::vector<std::string> bio_fasta_names(const std::string& path) {
    BioLineReader in(path);
    std::vector<std::string> out;
    std::string line;
    while (in.next(line)) {
        if (!line.empty() && line[0] == '>') {
            size_t a = 1;
            while (a < line.size() && std::isspace(static_cast<unsigned char>(line[a]))) ++a;
            size_t b = a;
            while (b < line.size() && !std::isspace(static_cast<unsigned char>(line[b]))) ++b;
            out.push_back(line.substr(a, b - a));
        }
    }
    return out;
}

// Write the spliced sequence of every transcript of `entries` as FASTA,
// reading the genome once; returns the transcript ids in the order written
// (genome order, then model order within a sequence). A transcript on a
// sequence absent from the genome is an error.
// Writes the spliced transcripts of the genome's sequences, as they come.
class BioTxWriter {
  public:
    BioTxWriter(const std::vector<BioGffEntry>& entries, const std::string& out) : w_(out) {
        std::unordered_set<std::string> seen;
        for (const auto& e : entries) {
            const std::string s = bio_entry_seqid(e);
            if (seen.insert(s).second) names_.push_back(s);
        }
        for (auto& m : bio_transcript_models(names_, entries)) bySeq_[m.seqid].push_back(std::move(m));
    }

    void add(const std::vector<BioFastaEntry>& batch) {
        for (const auto& g : batch) {
            const std::string name = bio_split(g.header + " ", ' ')[0];
            auto it = bySeq_.find(name);
            if (it == bySeq_.end()) continue;
            std::string text;
            for (const BioTxModel& m : it->second) {
                BioFastaEntry e{m.id, bio_extract(m.mapping, g.bases)};
                bio_fasta_text(text, e);
                written_.push_back(m.id);
            }
            w_.write(text);
            found_.insert(name);
        }
    }

    std::vector<std::string> finish(const std::string& genome) {
        for (const auto& n : names_) {
            if (!found_.count(n)) throw std::runtime_error("writeTranscriptFasta: sequence " + n + " of the annotation is not in " + genome);
        }
        w_.close();
        return written_;
    }

  private:
    BioWriter w_;
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<BioTxModel>> bySeq_;
    std::vector<std::string> written_;
    std::unordered_set<std::string> found_;
};

inline std::vector<std::string> bio_write_transcript_fasta(const std::string& genome, const std::vector<BioGffEntry>& entries,
                                                         const std::string& out) {
    BioTxWriter w(entries, out);
    bio_read_fasta(genome, [&](const std::vector<BioFastaEntry>& batch) { w.add(batch); });
    return w.finish(genome);
}

// Calls `f` on each batch of a pull until it is exhausted; a failed pull is
// an error.
template <class Pull, class F>
void bio_drain_pull(const std::string& what, Pull pull, F f) {
    for (;;) {
        auto r = pull();
        if (!std::get<0>(r)) throw std::runtime_error(what + ": " + std::get<1>(r));
        if (std::get<2>(r).empty()) return;
        f(std::get<2>(r));
    }
}

template <class Pull>
std::vector<std::string> bio_write_transcript_fasta_pull(const std::vector<BioGffEntry>& entries, const std::string& out, Pull pull) {
    BioTxWriter w(entries, out);
    bio_drain_pull("writeTranscriptFastaPull", pull, [&](const std::vector<BioFastaEntry>& batch) { w.add(batch); });
    return w.finish("the genome");
}

template <class Pull>
std::vector<std::string> bio_write_fasta_pull(const std::string& path, Pull pull) {
    BioWriter w(path);
    std::vector<std::string> names;
    bio_drain_pull("writeFastaPull", pull, [&](const std::vector<BioFastaEntry>& batch) {
        std::string text;
        for (const auto& e : batch) {
            names.push_back(bio_split(e.header + " ", ' ')[0]);
            bio_fasta_text(text, e);
        }
        w.write(text);
    });
    w.close();
    return names;
}

struct BioTxGene {
    std::string tx;
    std::string gene;
    std::string geneName;
    std::string biotype;
};

inline std::string bio_first_attr(const std::string& attrs, const std::vector<std::string>& keys) {
    const BioAttrs parsed = bio_attrs_parse(attrs);
    for (const auto& k : keys) {
        const std::vector<std::string>* v = bio_attr_get(parsed, k);
        if (v && !v->empty()) return (*v)[0];
    }
    return "";
}

// Each transcript's gene, gene name and biotype, in transcriptModels order
// of the trees' entries; a key is looked up on the gene, then the
// transcript, then the transcript's children.
inline std::vector<BioTxGene> bio_tx_genes(const std::vector<BioGffRoot>& roots) {
    const std::vector<BioGffEntry> entries = bio_gff_entries(roots);
    std::vector<std::string> names;
    std::unordered_set<std::string> seen;
    for (const auto& e : entries) {
        const std::string s = bio_entry_seqid(e);
        if (seen.insert(s).second) names.push_back(s);
    }
    auto hasExon = [](const BioFeature<BioGffNode>& f) {
        return std::any_of(f.children.begin(), f.children.end(), [](const BioFeature<BioGffNode>& c) { return c.node.kind == "exon"; });
    };
    std::unordered_map<std::string, std::vector<const std::string*>> search;
    for (const auto& r : roots) {
        const auto& g = r.feature;
        if (hasExon(g)) {
            std::vector<const std::string*> cols{&g.node.attrs};
            for (const auto& c : g.children) cols.push_back(&c.node.attrs);
            search[g.node.id] = cols;
        }
        for (const auto& t : g.children) {
            if (hasExon(t)) {
                std::vector<const std::string*> cols{&g.node.attrs, &t.node.attrs};
                for (const auto& c : t.children) cols.push_back(&c.node.attrs);
                search[t.node.id] = cols;
            }
        }
    }
    std::vector<BioTxGene> out;
    for (const auto& m : bio_transcript_models(names, entries)) {
        const auto it = search.find(m.id);
        auto pick = [&](const std::vector<std::string>& keys) -> std::string {
            if (it == search.end()) return "";
            for (const std::string* c : it->second) {
                std::string x = bio_first_attr(*c, keys);
                if (!x.empty()) return x;
            }
            return "";
        };
        const std::string name = pick({"gene_name", "Name"});
        out.push_back(BioTxGene{m.id, m.gene, name.empty() ? m.gene : name, pick({"gene_biotype", "gene_type", "biotype"})});
    }
    return out;
}


// Pulls over FASTQ text
// ---------------------
//
// A pull is a handle into a process-local table of open readers; `next`
// returns the next batch (empty at the end) and drops the reader then.

struct BioFastqPull {
    std::vector<std::vector<std::string>> runs;
    size_t next = 0;
    std::vector<std::unique_ptr<BioLineReader>> cur;
};

inline std::unordered_map<int64_t, BioFastqPull>& bio_pulls() {
    static std::unordered_map<int64_t, BioFastqPull> t;
    return t;
}

inline int64_t bio_open_pull(std::vector<std::vector<std::string>> runs) {
    static int64_t next = 1;
    const int64_t h = next++;
    bio_pulls()[h] = BioFastqPull{std::move(runs), 0, {}};
    return h;
}

inline int64_t bio_open_pairs(const std::vector<std::string>& r1, const std::vector<std::string>& r2) {
    if (r1.size() != r2.size()) {
        throw std::runtime_error("openPairs: " + std::to_string(r1.size()) + " first-mate files but " + std::to_string(r2.size()) + " second-mate files");
    }
    std::vector<std::vector<std::string>> runs;
    for (size_t i = 0; i < r1.size(); ++i) runs.push_back({r1[i], r2[i]});
    return bio_open_pull(std::move(runs));
}

inline int64_t bio_open_reads(const std::vector<std::string>& r) {
    std::vector<std::vector<std::string>> runs;
    for (const auto& x : r) runs.push_back({x});
    return bio_open_pull(std::move(runs));
}

// The next batch of entries, one per open file of the current run.
template <class Out>
std::tuple<bool, std::string, std::vector<Out>> bio_next_pull(int64_t h, size_t mates) {
    auto& table = bio_pulls();
    auto it = table.find(h);
    if (it == table.end()) return {true, "", {}};
    BioFastqPull& p = it->second;
    std::vector<Out> out;
    try {
        size_t weight = 0;
        while (weight < BIO_BATCH_BYTES) {
            if (p.cur.empty()) {
                if (p.next == p.runs.size()) break;
                for (const auto& f : p.runs[p.next]) p.cur.push_back(std::make_unique<BioLineReader>(f));
                ++p.next;
            }
            std::vector<BioFastqEntry> es(mates);
            size_t more = 0;
            for (size_t k = 0; k < mates; ++k) more += bio_next_fastq(*p.cur[k], es[k]);
            if (more == 0) {
                p.cur.clear();
                continue;
            }
            if (more != mates) {
                throw std::runtime_error("mate files " + p.runs[p.next - 1][0] + " and " + p.runs[p.next - 1][1] + " hold different numbers of reads");
            }
            for (const auto& e : es) weight += e.header.size() + 2 * e.bases.bases.size();
            if constexpr (std::is_same_v<Out, BioFastqEntry>) {
                out.push_back(std::move(es[0]));
            } else {
                out.push_back(Out{std::move(es[0]), std::move(es[1])});
            }
        }
    } catch (const std::exception& e) {
        table.erase(it);
        return {false, e.what(), {}};
    }
    if (out.empty()) table.erase(it);
    return {true, "", std::move(out)};
}

inline std::tuple<bool, std::string, std::vector<std::tuple<BioFastqEntry, BioFastqEntry>>> bio_next_pairs(int64_t h) {
    return bio_next_pull<std::tuple<BioFastqEntry, BioFastqEntry>>(h, 2);
}

inline std::tuple<bool, std::string, std::vector<BioFastqEntry>> bio_next_reads(int64_t h) {
    return bio_next_pull<BioFastqEntry>(h, 1);
}


// Mapper references on the genome
// -------------------------------

enum class BioLibStrand : uint8_t { LibUnstranded = 0, LibForward = 1, LibReverse = 2 };
enum class BioLayout : uint8_t { PairedEnd = 0, SingleEnd = 1 };

struct BioRefModels {
    std::vector<BioMapping<BioTranscript, BioGenome>> maps;
};

inline BioRefModels bio_ref_models(const std::vector<std::string>& names, const std::vector<BioTxModel>& models) {
    std::unordered_map<std::string, const BioTxModel*> byId;
    for (const auto& m : models) byId[m.id] = &m;
    BioRefModels out;
    for (const auto& n : names) {
        auto it = byId.find(n);
        if (it == byId.end()) {
            out.maps.push_back(BioMapping<BioTranscript, BioGenome>{-1, -1, {}});
        } else {
            out.maps.push_back(it->second->mapping);
        }
    }
    return out;
}

inline std::optional<std::tuple<int64_t, bool>> bio_project_point(const BioMapping<BioTranscript, BioGenome>& m, int64_t x) {
    if (x < 0) return std::nullopt;
    const uint64_t u = static_cast<uint64_t>(x);
    auto it = std::partition_point(m.blocks.begin(), m.blocks.end(), [&](const BioBlock& b) { return b.srcStart + b.len <= u; });
    if (it == m.blocks.end() || u < it->srcStart) return std::nullopt;
    const uint64_t off = u - it->srcStart;
    if (it->orient == BioOrient::Fwd) return std::make_tuple(static_cast<int64_t>(it->dstStart + off), true);
    return std::make_tuple(static_cast<int64_t>(it->dstStart + it->len - 1 - off), false);
}

#endif
