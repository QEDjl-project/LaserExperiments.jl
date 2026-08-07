struct LaserSystemMetadata
    status::LaserStatus
    institution::String
    location::String
    facility::String
    construction_date::Union{Date, Nothing}
    last_update::DateTime
end
