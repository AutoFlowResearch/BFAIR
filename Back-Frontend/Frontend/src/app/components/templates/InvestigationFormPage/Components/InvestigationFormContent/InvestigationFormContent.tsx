import React, { useReducer } from 'react';
import { Investigation } from '../../../../../shared/model/investigation.model';
import Card from '../../../../atoms/Card';
import InvestigationFormGeneral from '../InvestigationFormGeneral';
import InvestigationFormStudies from '../InvestigationFormStudies';
import backIcon from '../../../../../../assets/images/back-icon.svg';
import { Link } from 'react-router-dom';

const InvestigationFormContent = () => {
  return (
    <div className="investigation-form">
      <Card>
        <div className="investigation-form__inner">
          <h3 className="form-title">
            <Link to="/">
              <img src={backIcon} alt="" className="back-icon" />
            </Link>
            Investigation
          </h3>
          <InvestigationFormGeneral />
          <InvestigationFormStudies />
        </div>
        <div className="card__footer">
          <Link to="/">
            <button className="btn--flat">Cancel</button>
          </Link>
          <button className="btn--primary float--right">Submit</button>
        </div>
      </Card>
    </div>
  );
};

export default InvestigationFormContent;
