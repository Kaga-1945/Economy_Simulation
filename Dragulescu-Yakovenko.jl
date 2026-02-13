using Random
using Statistics
using StatsBase
using Plots
gr()

# エージェントの定義
mutable struct Agent
    # 識別(要らなかった)
    id::Int64
    # 所持金
    m::Float64
end

# エージェントの初期化
function init_agent(N::Int64, M::Float64)
    x = randexp(N)
    w = x ./ sum(x)
    agents = [Agent(i, M * w[i]) for i in 1:N]
    return agents
end

# 交換
function exchange(m1::Float64, m2::Float64)
    # 論文から
    s = m1 + m2
    ϵ = rand()
    return ϵ * s, (1 - ϵ) * s
end

# 取引
function deal!(agents::Vector{Agent}, steps::Int64)
    N = length(agents)
    # 規定の取引回数まで
    @inbounds @simd for _ in 1:steps
        # ランダムに2体選択
        i, j = rand(1:N, 2)
        # もし同じ2人を選んでしまった場合は再抽選する
        while i == j
            j = rand(1:N)
        end

        # 取引の実行
        m1, m2 = exchange(agents[i].m, agents[j].m)
        # 所持金の更新
        agents[i].m = m1
        agents[j].m = m2
    end
    return agents
end


function simlator(M::Float64, N::Int64, steps::Int64)
    # エージェントを初期化
    agents = init_agent(N, Float64(M))
    # 取引を実行
    deal!(agents, steps)
    return agents
end

function main()
    # お金の総額
    M = 1000.0
    # 総人数
    N = 1000
    # 取引回数
    steps = 2 * 10^7
    # シミュレーションの実行
    agents = simlator(M, N, steps)

    # 所持金をリスト化する
    m = [a.m for a in agents]
    # 所持金の平均を計算
    T = mean(m)

    m_sorted = sort(m)
    N = length(m_sorted)

    # 累積頻度分布を作る
    S = (N .- (1:N) .+ 1) ./ N

    # 理論：S(m) = exp(-m/T)
    mgrid = range(0, stop=maximum(m_sorted), length=400)
    S_theory = exp.(-mgrid ./ T)

    # 生存関数
    p1 = plot(m_sorted, S, seriestype=:scatter, markersize=2, markerstrokewidth=0,
        xlabel="m", ylabel="Pr(m_i ≥ m)", label="simulation",
        title="Survival function")
    plot!(p1, mgrid, S_theory, label="theory: exp(-m/T)")

    # 片対数：log S(m) vs m（指数なら直線）
    p2 = plot(m_sorted, log.(S), seriestype=:scatter, markersize=2, markerstrokewidth=0,
        xlabel="m", ylabel="log Pr(m_i ≥ m)", label="simulation",
        title="Semilog survival")
    plot!(p2, mgrid, -(mgrid ./ T), label="theory: -m/T")

    p = plot(p1, p2, layout=(1, 2), size=(1200, 450),
        eft_margin=10Plots.mm, bottom_margin=8Plots.mm, top_margin=6Plots.mm, right_margin=6Plots.mm)
    display(p)

end

main()