import React from 'react';

export const newField = ({
  input,
  placeholder,
  type,
  id,
  meta: { touched, error },
  ...rest
}) => {
  return (
    <div>
      <input
        {...input}
        type={type}
        placeholder={placeholder}
        id={id}
        className="form-control"
      />
      {touched && error && <p style={{ color: 'red' }}>{error}</p>}
    </div>
  );
};
